# src/scrna_benchmark/mixed_models.py

from pathlib import Path

import numpy as np
import pandas as pd

from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM

from .models import run_random_split_logreg, run_donor_split_logreg
from .splits import make_donor_folds


def logistic(x):
    """
    Logistic inverse-link function.
    """
    return 1.0 / (1.0 + np.exp(-x))


def result_to_prediction_metadata(
    result,
    obs_test=None,
    donor_col="patient_id",
    site_col=None,
    scheme=None,
    representation=None,
    rep_key=None,
    split_id=None,
    dataset_name=None,
):
    """
    Convert a model result dict into prediction-level metadata.

    Parameters
    ----------
    result : dict
        Output from run_random_split_logreg or run_donor_split_logreg.
        Must contain y_test and y_pred. Preferably contains obs_test.
    obs_test : pd.DataFrame or None
        Test-cell metadata. If None, uses result["obs_test"].
    donor_col : str
        Donor/sample column.
    site_col : str or None
        Optional site/batch/study column to retain.
    scheme : str or None
        Evaluation scheme label.
    representation : str or None
        Human-readable representation label.
    rep_key : str or None
        AnnData representation key.
    split_id : str or None
        Split identifier.
    dataset_name : str or None
        Dataset label.

    Returns
    -------
    pd.DataFrame
        Prediction-level metadata table.
    """
    if obs_test is None:
        if "obs_test" not in result:
            raise KeyError("obs_test was not provided and result does not contain 'obs_test'.")
        obs_test = result["obs_test"]

    if len(obs_test) != len(result["y_test"]):
        raise ValueError(
            "obs_test length does not match y_test length: "
            f"{len(obs_test)} != {len(result['y_test'])}"
        )

    cols = []

    if donor_col in obs_test.columns:
        cols.append(donor_col)
    else:
        raise KeyError(f"{donor_col} not found in obs_test.")

    if site_col is not None and site_col in obs_test.columns:
        cols.append(site_col)

    pred_df = obs_test[cols].copy()
    pred_df["y_true"] = result["y_test"]
    pred_df["y_pred"] = result["y_pred"]
    pred_df["correct"] = (pred_df["y_true"] == pred_df["y_pred"]).astype(int)

    if scheme is not None:
        pred_df["scheme"] = scheme

    if representation is not None:
        pred_df["representation"] = representation

    if rep_key is not None:
        pred_df["rep_key"] = rep_key

    if split_id is not None:
        pred_df["split_id"] = split_id

    if dataset_name is not None:
        pred_df["dataset"] = dataset_name

    return pred_df.reset_index(drop=True)


def collect_random_split_prediction_metadata(
    adata,
    representations,
    celltype_col="cell_type",
    donor_col="patient_id",
    site_col=None,
    batch_col=None,
    test_size=0.2,
    n_repeats=5,
    random_state=42,
    dataset_name=None,
    verbose=True,
):
    """
    Collect prediction-level metadata from repeated random cell-level splits.

    Returns
    -------
    pd.DataFrame
        One row per test-cell prediction.
    """
    all_pred_dfs = []

    for rep_label, rep_key in representations.items():
        if verbose:
            print(f"Collecting random split metadata for {rep_label} ({rep_key})")

        for repeat in range(n_repeats):
            seed = random_state + repeat

            res = run_random_split_logreg(
                adata=adata,
                rep_key=rep_key,
                celltype_col=celltype_col,
                batch_col=batch_col,
                test_size=test_size,
                random_state=seed,
            )

            pred_df = result_to_prediction_metadata(
                result=res,
                obs_test=res["obs_test"],
                donor_col=donor_col,
                site_col=site_col,
                scheme="random",
                representation=rep_label,
                rep_key=rep_key,
                split_id=f"random_{repeat + 1}",
                dataset_name=dataset_name,
            )

            all_pred_dfs.append(pred_df)

    if not all_pred_dfs:
        return pd.DataFrame()

    return pd.concat(all_pred_dfs, ignore_index=True)


def collect_donor_cv_prediction_metadata(
    adata,
    representations,
    celltype_col="cell_type",
    donor_col="patient_id",
    site_col=None,
    batch_col=None,
    n_folds=5,
    random_state=42,
    dataset_name=None,
    verbose=True,
):
    """
    Collect prediction-level metadata from donor-held-out CV.

    Returns
    -------
    pd.DataFrame
        One row per test-cell prediction.
    """
    donor_folds = make_donor_folds(
        adata=adata,
        donor_col=donor_col,
        n_folds=n_folds,
        random_state=random_state,
    )

    all_pred_dfs = []

    for rep_label, rep_key in representations.items():
        if verbose:
            print(f"Collecting donor-CV metadata for {rep_label} ({rep_key})")

        for fold_idx, test_donors in enumerate(donor_folds):
            train_donors = np.concatenate([
                donor_folds[j]
                for j in range(len(donor_folds))
                if j != fold_idx
            ])

            res = run_donor_split_logreg(
                adata=adata,
                rep_key=rep_key,
                train_donors=train_donors,
                test_donors=test_donors,
                celltype_col=celltype_col,
                donor_col=donor_col,
                batch_col=batch_col,
                random_state=random_state,
            )

            pred_df = result_to_prediction_metadata(
                result=res,
                obs_test=res["obs_test"],
                donor_col=donor_col,
                site_col=site_col,
                scheme="donor_held_out",
                representation=rep_label,
                rep_key=rep_key,
                split_id=f"fold_{fold_idx + 1}",
                dataset_name=dataset_name,
            )

            all_pred_dfs.append(pred_df)

    if not all_pred_dfs:
        return pd.DataFrame()

    return pd.concat(all_pred_dfs, ignore_index=True)


def collect_scheme_prediction_metadata(
    adata,
    representations,
    celltype_col="cell_type",
    donor_col="patient_id",
    site_col=None,
    batch_col=None,
    test_size=0.2,
    n_random_repeats=5,
    n_folds=5,
    random_state=42,
    dataset_name=None,
    verbose=True,
):
    """
    Collect prediction-level metadata for both:
    1. repeated random cell-level splits
    2. donor-held-out CV

    Returns
    -------
    pd.DataFrame
        Combined prediction-level metadata.
    """
    random_df = collect_random_split_prediction_metadata(
        adata=adata,
        representations=representations,
        celltype_col=celltype_col,
        donor_col=donor_col,
        site_col=site_col,
        batch_col=batch_col,
        test_size=test_size,
        n_repeats=n_random_repeats,
        random_state=random_state,
        dataset_name=dataset_name,
        verbose=verbose,
    )

    donor_df = collect_donor_cv_prediction_metadata(
        adata=adata,
        representations=representations,
        celltype_col=celltype_col,
        donor_col=donor_col,
        site_col=site_col,
        batch_col=batch_col,
        n_folds=n_folds,
        random_state=random_state,
        dataset_name=dataset_name,
        verbose=verbose,
    )

    return pd.concat([random_df, donor_df], ignore_index=True)


def fit_correctness_glmm(
    pred_df,
    donor_col="patient_id",
    scheme_col="scheme",
    outcome_col="correct",
    reference_scheme="donor_held_out",
):
    """
    Fit a binomial mixed model:

        correct ~ C(scheme) + donor random intercept

    Parameters
    ----------
    pred_df : pd.DataFrame
        Prediction-level metadata table.
    donor_col : str
        Donor/sample column.
    scheme_col : str
        Column containing scheme labels.
    outcome_col : str
        Binary correctness column.
    reference_scheme : str
        Reference category for scheme.

    Returns
    -------
    model, result
    """
    df = pred_df.copy()

    required_cols = [donor_col, scheme_col, outcome_col]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns for GLMM: {missing}")

    df[donor_col] = df[donor_col].astype("category")

    observed_schemes = sorted(df[scheme_col].astype(str).unique())

    if reference_scheme not in observed_schemes:
        raise ValueError(
            f"reference_scheme={reference_scheme} not found. "
            f"Observed schemes: {observed_schemes}"
        )

    scheme_order = [reference_scheme] + [
        scheme for scheme in observed_schemes
        if scheme != reference_scheme
    ]

    df[scheme_col] = pd.Categorical(
        df[scheme_col].astype(str),
        categories=scheme_order,
        ordered=False,
    )

    formula = f"{outcome_col} ~ C({scheme_col})"
    vc_formulas = {
        "donor_re": f"0 + C({donor_col})"
    }

    model = BinomialBayesMixedGLM.from_formula(
        formula,
        vc_formulas,
        df,
    )

    result = model.fit_vb()

    return model, result


def summarize_correctness_glmm(
    model,
    result,
    representation=None,
    reference_scheme="donor_held_out",
    comparison_scheme="random",
):
    """
    Summarize GLMM fixed effects as predicted probabilities.

    Returns
    -------
    summary_df : pd.DataFrame
        Predicted probability correct for reference and comparison schemes.
    coef_df : pd.DataFrame
        Fixed-effect coefficient table.
    """
    fe_names = model.exog_names
    fe_mean = pd.Series(result.fe_mean, index=fe_names)
    fe_sd = pd.Series(result.fe_sd, index=fe_names)

    if "Intercept" not in fe_mean.index:
        raise KeyError("GLMM fixed effects do not contain an Intercept.")

    intercept = fe_mean["Intercept"]
    p_reference = logistic(intercept)

    scheme_terms = [
        name for name in fe_names
        if "C(" in name and comparison_scheme in name
    ]

    if len(scheme_terms) == 0:
        raise KeyError(
            f"No fixed-effect term found for comparison_scheme={comparison_scheme}. "
            f"Available fixed effects: {fe_names}"
        )

    scheme_coef_name = scheme_terms[0]
    scheme_coef = fe_mean[scheme_coef_name]
    p_comparison = logistic(intercept + scheme_coef)

    summary_df = pd.DataFrame([
        {
            "representation": representation,
            "scheme": reference_scheme,
            "pred_prob_correct": p_reference,
            "coef_intercept": intercept,
            "coef_scheme_comparison": scheme_coef,
            "scheme_coef_name": scheme_coef_name,
        },
        {
            "representation": representation,
            "scheme": comparison_scheme,
            "pred_prob_correct": p_comparison,
            "coef_intercept": intercept,
            "coef_scheme_comparison": scheme_coef,
            "scheme_coef_name": scheme_coef_name,
        },
    ])

    coef_df = pd.DataFrame({
        "term": fe_names,
        "mean": fe_mean.values,
        "sd": fe_sd.values,
    })

    return summary_df, coef_df


def fit_glmm_by_representation(
    pred_df,
    donor_col="patient_id",
    representation_col="representation",
    scheme_col="scheme",
    outcome_col="correct",
    reference_scheme="donor_held_out",
    comparison_scheme="random",
    verbose=True,
):
    """
    Fit separate correctness GLMMs for each representation.

    Returns
    -------
    summary_all : pd.DataFrame
        Predicted probabilities by representation and scheme.
    coef_parts : dict[str, pd.DataFrame]
        Fixed-effect coefficient tables by representation.
    fitted : dict
        Raw fitted model/result objects by representation.
    """
    if representation_col not in pred_df.columns:
        raise KeyError(f"{representation_col} not found in pred_df.")

    summary_parts = []
    coef_parts = {}
    fitted = {}

    reps = sorted(pred_df[representation_col].astype(str).unique())

    for rep in reps:
        if verbose:
            print(f"Fitting correctness GLMM for {rep}")

        df_rep = pred_df.loc[
            pred_df[representation_col].astype(str) == rep
        ].copy()

        model, result = fit_correctness_glmm(
            pred_df=df_rep,
            donor_col=donor_col,
            scheme_col=scheme_col,
            outcome_col=outcome_col,
            reference_scheme=reference_scheme,
        )

        summary_df, coef_df = summarize_correctness_glmm(
            model=model,
            result=result,
            representation=rep,
            reference_scheme=reference_scheme,
            comparison_scheme=comparison_scheme,
        )

        summary_parts.append(summary_df)
        coef_parts[rep] = coef_df
        fitted[rep] = {
            "model": model,
            "result": result,
        }

    summary_all = pd.concat(summary_parts, ignore_index=True)

    return summary_all, coef_parts, fitted


def save_glmm_outputs(
    pred_df,
    summary_df,
    coef_parts,
    out_dir,
    prefix="glmm",
):
    """
    Save GLMM input predictions and summaries.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pred_df.to_csv(out_dir / f"{prefix}_input_predictions.csv", index=False)
    summary_df.to_csv(out_dir / f"{prefix}_predicted_probabilities.csv", index=False)

    for rep, coef_df in coef_parts.items():
        safe_rep = str(rep).replace("/", "_").replace(" ", "_")
        coef_df.to_csv(
            out_dir / f"{prefix}_{safe_rep}_fixed_effects.csv",
            index=False,
        )