# src/scrna_benchmark/plotting.py

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.colors import PowerNorm
from sklearn.metrics import confusion_matrix


DEFAULT_REP_ORDER = ["hvg", "pca", "harmony", "scvi"]

DEFAULT_REP_LABELS = {
    "hvg": "HVG",
    "pca": "PCA",
    "harmony": "Harmony",
    "scvi": "scVI",
}

REP_ALIASES = {
    "hvg": "hvg",
    "pca": "pca",
    "X_pca": "pca",
    "harmony": "harmony",
    "X_harmony": "harmony",
    "X_pca_harmony": "harmony",
    "scvi": "scvi",
    "X_scVI": "scvi",
}


def canonicalize_representation(x):
    """Map representation names/keys to a standard plotting label."""
    return REP_ALIASES.get(str(x), str(x))


def _read_df(x):
    """Accept either a DataFrame or CSV path."""
    if isinstance(x, (str, Path)):
        return pd.read_csv(x)
    return x.copy()


def _save_figure(fig, out_file):
    if out_file is None:
        return

    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=300, bbox_inches="tight")


def _annotate_bars(ax, bars, offset=0.005, fmt=".3f"):
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h + offset,
            format(h, fmt),
            ha="center",
            va="bottom",
            fontsize=9,
        )


def _summarize_repeated_metrics(df):
    """
    Summarize repeated random-split results.

    Expects one row per representation x repeat.
    """
    df = _read_df(df)
    df["representation"] = df["representation"].map(
        canonicalize_representation
    )

    return (
        df
        .groupby("representation", as_index=False)
        .agg(
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
            accuracy_mean=("accuracy", "mean"),
            accuracy_std=("accuracy", "std"),
        )
    )

# main figure
def plot_random_vs_donor_cv(
    random_metrics,
    donor_metrics,
    rep_order=DEFAULT_REP_ORDER,
    rep_labels=DEFAULT_REP_LABELS,
    out_file=None,
    title="Random split overestimates cross-donor performance",
):
    """
    Compare repeated random cell-level splits with donor-held-out CV
    using macro F1.

    Parameters
    ----------
    random_metrics
        DataFrame or path to random_split_repeated_metrics.csv.

    donor_metrics
        DataFrame or path to donor_cv/metrics.csv.

    rep_order
        Canonical representation order.

    rep_labels
        Mapping from canonical representation name to display label.

    out_file
        Optional figure output path.

    title
        Plot title.
    """
    random_df = _summarize_repeated_metrics(random_metrics)

    donor_df = _read_df(donor_metrics)
    donor_df["representation"] = donor_df["representation"].map(
        canonicalize_representation
    )

    common_reps = [
        rep for rep in rep_order
        if (
            rep in set(random_df["representation"])
            and rep in set(donor_df["representation"])
        )
    ]

    random_plot = (
        random_df
        .set_index("representation")
        .reindex(common_reps)
    )

    donor_plot = (
        donor_df
        .set_index("representation")
        .reindex(common_reps)
    )

    x = np.arange(len(common_reps))
    width = 0.36

    fig, ax = plt.subplots(figsize=(8, 5))

    bars_random = ax.bar(
        x - width / 2,
        random_plot["macro_f1_mean"],
        width,
        yerr=random_plot["macro_f1_std"],
        capsize=5,
        label="Random cell-level split",
        edgecolor="black",
    )

    bars_donor = ax.bar(
        x + width / 2,
        donor_plot["macro_f1_mean"],
        width,
        yerr=donor_plot["macro_f1_std"],
        capsize=5,
        label="Donor-held-out CV",
        edgecolor="black",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [rep_labels[r] for r in common_reps]
    )

    ax.set_xlabel("Representation")
    ax.set_ylabel("Macro F1")
    ax.set_title(title)
    ax.legend(frameon=False)

    lower = min(
        (
            random_plot["macro_f1_mean"]
            - random_plot["macro_f1_std"]
        ).min(),
        (
            donor_plot["macro_f1_mean"]
            - donor_plot["macro_f1_std"]
        ).min(),
    )

    upper = max(
        (
            random_plot["macro_f1_mean"]
            + random_plot["macro_f1_std"]
        ).max(),
        (
            donor_plot["macro_f1_mean"]
            + donor_plot["macro_f1_std"]
        ).max(),
    )

    ax.set_ylim(
        max(0, lower - 0.03),
        min(1, upper + 0.04),
    )

    _annotate_bars(ax, bars_random, offset=0.006)
    _annotate_bars(ax, bars_donor, offset=0.006)

    fig.tight_layout()
    _save_figure(fig, out_file)

    return fig, ax


# ============================================================
# Batch-covariate sensitivity
# ============================================================

def plot_batch_covariate_comparison(
    no_batch_metrics,
    with_batch_metrics,
    rep_order=DEFAULT_REP_ORDER,
    rep_labels=DEFAULT_REP_LABELS,
    out_file=None,
    title="Effect of batch covariate on donor-held-out performance",
):
    """
    Compare donor-held-out macro F1 with and without a batch covariate.
    """
    no_batch = _read_df(no_batch_metrics)
    with_batch = _read_df(with_batch_metrics)

    for df in (no_batch, with_batch):
        df["representation"] = df["representation"].map(
            canonicalize_representation
        )

    no_batch = no_batch.set_index("representation").reindex(rep_order)
    with_batch = with_batch.set_index("representation").reindex(rep_order)

    x = np.arange(len(rep_order))
    width = 0.36

    fig, ax = plt.subplots(figsize=(8, 5))

    bars_no = ax.bar(
        x - width / 2,
        no_batch["macro_f1_mean"],
        width,
        yerr=no_batch["macro_f1_std"],
        capsize=5,
        label="No batch covariate",
        edgecolor="black",
    )

    bars_yes = ax.bar(
        x + width / 2,
        with_batch["macro_f1_mean"],
        width,
        yerr=with_batch["macro_f1_std"],
        capsize=5,
        label="With batch covariate",
        edgecolor="black",
    )

    ax.set_xticks(x)
    ax.set_xticklabels([rep_labels[r] for r in rep_order])
    ax.set_xlabel("Representation")
    ax.set_ylabel("Macro F1")
    ax.set_title(title)
    ax.legend(frameon=False)

    lower = min(
        (no_batch["macro_f1_mean"] - no_batch["macro_f1_std"]).min(),
        (with_batch["macro_f1_mean"] - with_batch["macro_f1_std"]).min(),
    )
    upper = max(
        (no_batch["macro_f1_mean"] + no_batch["macro_f1_std"]).max(),
        (with_batch["macro_f1_mean"] + with_batch["macro_f1_std"]).max(),
    )

    ax.set_ylim(
        max(0, lower - 0.03),
        min(1, upper + 0.04),
    )

    _annotate_bars(ax, bars_no)
    _annotate_bars(ax, bars_yes)

    fig.tight_layout()
    _save_figure(fig, out_file)

    return fig, ax


def _coerce_macro_f1_summary(metrics):
    """Return representation-level macro-F1 mean/std from raw or summarized metrics."""
    df = _read_df(metrics)
    if "representation" not in df.columns:
        raise ValueError("Metrics table must contain a 'representation' column.")
    df["representation"] = df["representation"].map(canonicalize_representation)

    if {"macro_f1_mean", "macro_f1_std"}.issubset(df.columns):
        out = df[["representation", "macro_f1_mean", "macro_f1_std"]].copy()
        return out.drop_duplicates("representation")

    if "macro_f1" not in df.columns:
        raise ValueError(
            "Metrics table must contain either macro_f1 or "
            "macro_f1_mean/macro_f1_std."
        )

    return (
        df.groupby("representation", as_index=False)
        .agg(
            macro_f1_mean=("macro_f1", "mean"),
            macro_f1_std=("macro_f1", "std"),
        )
        .fillna({"macro_f1_std": 0.0})
    )


def _infer_random_metrics_from_donor_path(donor_metrics):
    """Infer a sibling random-split metrics file from a donor-CV metrics path."""
    if not isinstance(donor_metrics, (str, Path)):
        raise ValueError(
            "Automatic random-split path inference requires donor metrics to be "
            "passed as a CSV path. If using DataFrames, pass the random metrics "
            "explicitly via random_no_batch_metrics/random_with_batch_metrics."
        )

    donor_path = Path(donor_metrics)
    if donor_path.parent.name != "donor_cv":
        raise ValueError(
            f"Could not infer result root from {donor_path}. Expected a path "
            "inside a donor_cv/ directory."
        )

    random_dir = donor_path.parent.parent / "random_split"
    candidates = [
        random_dir / "random_split_repeated_metrics.csv",
        random_dir / "metrics.csv",
    ]
    candidates.extend(sorted(random_dir.glob("*repeated*metrics*.csv")))

    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "Could not find random-split metrics next to donor-CV results. Looked in "
        f"{random_dir}. Pass the file explicitly if it has a nonstandard name."
    )


def plot_batch_covariate_by_scheme(
    no_batch_metrics,
    with_batch_metrics,
    random_no_batch_metrics=None,
    random_with_batch_metrics=None,
    rep_order=DEFAULT_REP_ORDER,
    rep_labels=DEFAULT_REP_LABELS,
    out_file=None,
    title="Batch/source covariate sensitivity across evaluation schemes",
):
    """Compare no-covariate vs covariate models *within both evaluation schemes*.

    Parameters
    ----------
    no_batch_metrics, with_batch_metrics
        Donor-CV summary tables or CSV paths. These are intentionally the first
        two arguments so existing figure cells can be migrated with a minimal
        function-name change.
    random_no_batch_metrics, random_with_batch_metrics
        Repeated random-split metric tables or paths. If omitted and the donor
        metrics are paths inside ``<result_root>/donor_cv/``, the matching
        ``<result_root>/random_split/`` files are inferred automatically.

    Notes
    -----
    The shared y-axis and side-by-side evaluation panels are deliberate: the
    visual question is whether adding the batch/source covariate changes the
    random-vs-donor evaluation story, not merely whether donor-CV F1 moves.
    """
    if random_no_batch_metrics is None:
        random_no_batch_metrics = _infer_random_metrics_from_donor_path(
            no_batch_metrics
        )
    if random_with_batch_metrics is None:
        random_with_batch_metrics = _infer_random_metrics_from_donor_path(
            with_batch_metrics
        )

    summaries = {
        "random_no": _coerce_macro_f1_summary(random_no_batch_metrics),
        "random_yes": _coerce_macro_f1_summary(random_with_batch_metrics),
        "donor_no": _coerce_macro_f1_summary(no_batch_metrics),
        "donor_yes": _coerce_macro_f1_summary(with_batch_metrics),
    }

    rep_sets = [set(df["representation"]) for df in summaries.values()]
    common_reps = [rep for rep in rep_order if all(rep in s for s in rep_sets)]
    if not common_reps:
        raise ValueError("No common representations are present in all four inputs.")

    indexed = {
        key: df.set_index("representation").reindex(common_reps)
        for key, df in summaries.items()
    }

    all_lows = []
    all_highs = []
    for df in indexed.values():
        mean = df["macro_f1_mean"].to_numpy(dtype=float)
        std = df["macro_f1_std"].fillna(0).to_numpy(dtype=float)
        all_lows.extend((mean - std).tolist())
        all_highs.extend((mean + std).tolist())
    lower = max(0.0, float(np.nanmin(all_lows)) - 0.035)
    upper = min(1.0, float(np.nanmax(all_highs)) + 0.075)

    x = np.arange(len(common_reps))
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), sharey=True)

    panel_specs = [
        (axes[0], "random_no", "random_yes", "A. Random cell-level split"),
        (axes[1], "donor_no", "donor_yes", "B. Donor-held-out CV"),
    ]

    legend_handles = None
    legend_labels = None
    for ax, key_no, key_yes, panel_title in panel_specs:
        no_df = indexed[key_no]
        yes_df = indexed[key_yes]
        bars_no = ax.bar(
            x - width / 2,
            no_df["macro_f1_mean"],
            width,
            yerr=no_df["macro_f1_std"].fillna(0),
            capsize=4,
            label="No batch/source covariate",
            edgecolor="black",
        )
        bars_yes = ax.bar(
            x + width / 2,
            yes_df["macro_f1_mean"],
            width,
            yerr=yes_df["macro_f1_std"].fillna(0),
            capsize=4,
            label="With batch/source covariate",
            edgecolor="black",
        )
        ax.set_xticks(x)
        ax.set_xticklabels([rep_labels.get(r, r) for r in common_reps])
        ax.set_xlabel("Representation")
        ax.set_title(panel_title)
        ax.set_ylim(lower, upper)

        # Annotate the quantity the reader actually cares about here: how much
        # adding the covariate changes F1 within the same evaluation strategy.
        for i, rep in enumerate(common_reps):
            no_mean = float(no_df.loc[rep, "macro_f1_mean"])
            yes_mean = float(yes_df.loc[rep, "macro_f1_mean"])
            no_sd_value = no_df.loc[rep, "macro_f1_std"]
            yes_sd_value = yes_df.loc[rep, "macro_f1_std"]
            no_sd = 0.0 if pd.isna(no_sd_value) else float(no_sd_value)
            yes_sd = 0.0 if pd.isna(yes_sd_value) else float(yes_sd_value)
            y = max(no_mean + no_sd, yes_mean + yes_sd) + 0.012
            delta = yes_mean - no_mean
            ax.text(i, y, f"Δ {delta:+.3f}", ha="center", va="bottom", fontsize=8)

        if legend_handles is None:
            legend_handles = [bars_no, bars_yes]
            legend_labels = [
                "No batch/source covariate",
                "With batch/source covariate",
            ]

    axes[0].set_ylabel("Macro F1")
    fig.suptitle(title, y=1.02)
    fig.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        frameon=False,
        ncol=2,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    _save_figure(fig, out_file)
    return fig, axes

# ============================================================
# Donor ablation
# ============================================================

def plot_donor_ablation(
    summary_df,
    rep_order=DEFAULT_REP_ORDER,
    rep_labels=DEFAULT_REP_LABELS,
    out_file=None,
    title="Donor ablation under donor-held-out evaluation",
):
    """Plot macro F1 versus number of training donors."""
    df = _read_df(summary_df)
    df["representation"] = df["representation"].map(
        canonicalize_representation
    )

    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    for rep in rep_order:
        sub = (
            df[df["representation"] == rep]
            .sort_values("k_train_donors")
            .copy()
        )

        ax.plot(
            sub["k_train_donors"],
            sub["macro_f1_mean"],
            marker="o",
            linewidth=2,
            label=rep_labels[rep],
        )

        ax.errorbar(
            sub["k_train_donors"],
            sub["macro_f1_mean"],
            yerr=sub["macro_f1_std"],
            fmt="none",
            capsize=4,
        )

    ax.set_xlabel("Number of training donors")
    ax.set_ylabel("Macro F1")
    ax.set_title(title)
    ax.legend(frameon=False, title="Representation")

    ax.set_xticks(sorted(df["k_train_donors"].unique()))

    lower = (df["macro_f1_mean"] - df["macro_f1_std"]).min()
    upper = (df["macro_f1_mean"] + df["macro_f1_std"]).max()

    ax.set_ylim(
        max(0, lower - 0.03),
        min(1, upper + 0.03),
    )

    fig.tight_layout()
    _save_figure(fig, out_file)

    return fig, ax


# ============================================================
# Per-class F1
# ============================================================

def load_per_class_f1_matrix(
    donor_cv_dir,
    rep_order=DEFAULT_REP_ORDER,
    scheme_label="donor_cv",
):
    """
    Load mean per-class F1 tables saved by run_donor_cv_experiment().

    Representations without an output file are skipped. This allows datasets
    to use different valid representation sets (e.g. Blood Atlas has no scVI).
    """
    donor_cv_dir = Path(donor_cv_dir)

    parts = []
    available_reps = []

    for rep in rep_order:
        path = (
            donor_cv_dir
            / f"{scheme_label}_{rep}_mean_per_class_f1.csv"
        )

        if not path.exists():
            continue

        df = pd.read_csv(path)

        series = (
            df
            .set_index("cell_type")["f1"]
            .rename(rep)
        )

        parts.append(series)
        available_reps.append(rep)

    if not parts:
        raise FileNotFoundError(
            f"No per-class F1 files found in {donor_cv_dir} "
            f"for representations: {list(rep_order)}"
        )

    matrix = pd.concat(parts, axis=1)

    matrix["mean_f1"] = (
        matrix[available_reps]
        .mean(axis=1)
    )

    matrix = (
        matrix
        .sort_values(
            "mean_f1",
            ascending=False,
        )
        .drop(columns="mean_f1")
    )

    return matrix


def plot_per_class_f1_heatmap(
    f1_matrix,
    rep_order=DEFAULT_REP_ORDER,
    rep_labels=DEFAULT_REP_LABELS,
    annotate=True,
    out_file=None,
    title="Per-class F1 under donor-held-out evaluation",
):
    """Plot cell type x representation mean F1 heatmap."""
    matrix = f1_matrix.copy()

    # Keep only representations actually present in this dataset.
    available_reps = [
        rep
        for rep in rep_order
        if rep in matrix.columns
    ]

    if not available_reps:
        raise ValueError(
            "No requested representations are present in f1_matrix."
        )

    fig_height = max(
        6,
        0.45 * len(matrix),
    )

    fig, ax = plt.subplots(
        figsize=(7.5, fig_height)
    )

    im = ax.imshow(
        matrix[available_reps].values,
        aspect="auto",
        vmin=0,
        vmax=1,
        cmap="viridis",
    )

    ax.set_xticks(
        np.arange(
            len(available_reps)
        )
    )

    ax.set_xticklabels(
        [
            rep_labels[r]
            for r in available_reps
        ]
    )

    ax.set_yticks(
        np.arange(
            len(matrix.index)
        )
    )

    ax.set_yticklabels(
        matrix.index
    )

    ax.set_xlabel(
        "Representation"
    )
    ax.set_ylabel(
        "Cell type"
    )
    ax.set_title(
        title
    )

    cbar = fig.colorbar(
        im,
        ax=ax,
    )
    cbar.set_label(
        "Mean per-class F1"
    )

    if annotate:
        for i in range(
            matrix.shape[0]
        ):
            for j, rep in enumerate(
                available_reps
            ):
                val = matrix.iloc[i][rep]

                ax.text(
                    j,
                    i,
                    f"{val:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=(
                        "white"
                        if val < 0.45
                        else "black"
                    ),
                )

    fig.tight_layout()
    _save_figure(
        fig,
        out_file,
    )

    return fig, ax

# ============================================================
# Pooled confusion matrix
# ============================================================

def plot_confusion_matrix_from_predictions(
    predictions,
    out_file=None,
    title="Confusion matrix",
    gamma=0.5,
):
    """
    Build a confusion matrix from pooled prediction rows.

    Particularly useful for donor-CV all_predictions.csv.
    """
    pred = _read_df(predictions)

    labels = sorted(
        set(pred["y_true"].astype(str))
        | set(pred["y_pred"].astype(str))
    )

    cm = confusion_matrix(
        pred["y_true"].astype(str),
        pred["y_pred"].astype(str),
        labels=labels,
    )

    row_sums = cm.sum(axis=1, keepdims=True)

    cm_norm = np.divide(
        cm.astype(float),
        row_sums,
        out=np.zeros_like(cm, dtype=float),
        where=row_sums != 0,
    )

    fig, ax = plt.subplots(figsize=(10, 9))

    im = ax.imshow(
        cm_norm,
        cmap="viridis",
        norm=PowerNorm(
            gamma=gamma,
            vmin=0,
            vmax=1,
        ),
    )

    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))

    ax.set_xticklabels(labels, rotation=90)
    ax.set_yticklabels(labels)

    ax.set_xlabel("Predicted label")
    ax.set_ylabel("True label")
    ax.set_title(title)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Row-normalized proportion")

    fig.tight_layout()
    _save_figure(fig, out_file)

    cm_df = pd.DataFrame(
        cm_norm,
        index=labels,
        columns=labels,
    )

    return fig, ax, cm_df


# ============================================================
# Cross-group transfer
# ============================================================

def plot_group_transfer(
    summary_df,
    rep_order=DEFAULT_REP_ORDER,
    rep_labels=DEFAULT_REP_LABELS,
    group_label_map=None,
    out_file=None,
    title="Cross-group generalization",
):
    """
    Plot directional group/site transfer performance.
    """
    df = _read_df(summary_df)

    df["representation"] = df["representation"].map(
        canonicalize_representation
    )

    # Protect against accidentally duplicated bidirectional runs.
    if {"train_group", "test_group"}.issubset(df.columns):
        df = df.drop_duplicates(
            subset=[
                "representation",
                "train_group",
                "test_group",
            ]
        )

        label_map = group_label_map or {}

        def direction_label(row):
            train = label_map.get(
                str(row["train_group"]),
                str(row["train_group"]),
            )
            test = label_map.get(
                str(row["test_group"]),
                str(row["test_group"]),
            )
            return f"{train} → {test}"

        df["direction_label"] = df.apply(
            direction_label,
            axis=1,
        )

    else:
        df["direction_label"] = df["direction"]

    directions = list(
        dict.fromkeys(df["direction_label"])
    )

    pivot = (
        df
        .pivot(
            index="representation",
            columns="direction_label",
            values="macro_f1",
        )
        .reindex(rep_order)
    )

    x = np.arange(len(rep_order))
    n_directions = len(directions)
    width = 0.8 / n_directions

    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    all_bars = []

    for i, direction in enumerate(directions):
        offset = (
            i - (n_directions - 1) / 2
        ) * width

        bars = ax.bar(
            x + offset,
            pivot[direction],
            width,
            label=direction,
            edgecolor="black",
        )

        all_bars.append(bars)

    ax.set_xticks(x)
    ax.set_xticklabels([rep_labels[r] for r in rep_order])

    ax.set_xlabel("Representation")
    ax.set_ylabel("Macro F1")
    ax.set_title(title)
    ax.legend(frameon=False)

    lower = pivot.min().min()
    upper = pivot.max().max()

    ax.set_ylim(
        max(0, lower - 0.04),
        min(1, upper + 0.05),
    )

    for bars in all_bars:
        _annotate_bars(ax, bars)

    fig.tight_layout()
    _save_figure(fig, out_file)

    return fig, ax


# ============================================================
# GLMM sensitivity figure
# ============================================================

def plot_glmm_scheme_comparison(
    random_metrics,
    donor_metrics,
    glmm_summary,
    rep_order=DEFAULT_REP_ORDER,
    rep_labels=DEFAULT_REP_LABELS,
    out_file=None,
    title="Robustness across accuracy and mixed-effects modeling",
):
    """
    Two-panel figure:
      A. observed accuracy, random vs donor-held-out
      B. GLMM estimated probability correct
    """
    random_df = _summarize_repeated_metrics(random_metrics)

    donor_df = _read_df(donor_metrics)
    donor_df["representation"] = donor_df["representation"].map(
        canonicalize_representation
    )

    glmm_df = _read_df(glmm_summary)
    glmm_df["representation"] = glmm_df["representation"].map(
        canonicalize_representation
    )

    random_plot = (
        random_df
        .set_index("representation")
        .reindex(rep_order)
    )

    donor_plot = (
        donor_df
        .set_index("representation")
        .reindex(rep_order)
    )

    random_glmm = (
        glmm_df[glmm_df["scheme"] == "random"]
        .set_index("representation")
        .reindex(rep_order)
    )

    donor_glmm = (
        glmm_df[
            glmm_df["scheme"] == "donor_held_out"
        ]
        .set_index("representation")
        .reindex(rep_order)
    )

    x = np.arange(len(rep_order))
    width = 0.36

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13, 5.2),
    )

    # ------------------------
    # Panel A: raw accuracy
    # ------------------------
    bars_random = axes[0].bar(
        x - width / 2,
        random_plot["accuracy_mean"],
        width,
        yerr=random_plot["accuracy_std"],
        capsize=5,
        label="Random split",
        edgecolor="black",
    )

    bars_donor = axes[0].bar(
        x + width / 2,
        donor_plot["accuracy_mean"],
        width,
        yerr=donor_plot["accuracy_std"],
        capsize=5,
        label="Donor-held-out",
        edgecolor="black",
    )

    axes[0].set_xticks(x)
    axes[0].set_xticklabels(
        [rep_labels[r] for r in rep_order]
    )
    axes[0].set_ylabel("Accuracy")
    axes[0].set_title("A. Accuracy across evaluation schemes")
    axes[0].legend(frameon=False)

    acc_low = min(
        (
            random_plot["accuracy_mean"]
            - random_plot["accuracy_std"]
        ).min(),
        (
            donor_plot["accuracy_mean"]
            - donor_plot["accuracy_std"]
        ).min(),
    )

    acc_high = max(
        (
            random_plot["accuracy_mean"]
            + random_plot["accuracy_std"]
        ).max(),
        (
            donor_plot["accuracy_mean"]
            + donor_plot["accuracy_std"]
        ).max(),
    )

    axes[0].set_ylim(
        max(0, acc_low - 0.03),
        min(1, acc_high + 0.04),
    )

    _annotate_bars(axes[0], bars_random, offset=0.004)
    _annotate_bars(axes[0], bars_donor, offset=0.004)

    # ------------------------
    # Panel B: GLMM
    # ------------------------
    bars_random_glmm = axes[1].bar(
        x - width / 2,
        random_glmm["pred_prob_correct"],
        width,
        label="Random split",
        edgecolor="black",
    )

    bars_donor_glmm = axes[1].bar(
        x + width / 2,
        donor_glmm["pred_prob_correct"],
        width,
        label="Donor-held-out",
        edgecolor="black",
    )

    axes[1].set_xticks(x)
    axes[1].set_xticklabels(
        [rep_labels[r] for r in rep_order]
    )

    axes[1].set_ylabel(
        "Estimated probability of correct classification"
    )
    axes[1].set_title(
        "B. GLMM estimated probability"
    )
    axes[1].legend(frameon=False)

    glmm_low = glmm_df["pred_prob_correct"].min()
    glmm_high = glmm_df["pred_prob_correct"].max()

    axes[1].set_ylim(
        max(0, glmm_low - 0.03),
        min(1, glmm_high + 0.04),
    )

    _annotate_bars(
        axes[1],
        bars_random_glmm,
        offset=0.004,
    )
    _annotate_bars(
        axes[1],
        bars_donor_glmm,
        offset=0.004,
    )

    fig.suptitle(title, y=1.02)

    fig.tight_layout()
    _save_figure(fig, out_file)

    return fig, axes