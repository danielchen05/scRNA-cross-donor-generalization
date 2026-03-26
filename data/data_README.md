This project uses datasets accessed through the `pertpy` package through the commands:  
```python
import pertpy as pt
adata = pt.data.stephenson_2021_subsampled()
# adata.write("../data/stephenson_2021_subsampled.h5ad")
```
Data could be read in in two ways. We will provide code to read in data in both ways in each notebook/script, with default being option 1 (using `pertpy`).   
(1) Using the commands above.   
(2) Raw data is also saved as a `AnnData` object to avoid dependency of `pertpy` package (it works weird on Windows laptops). The `AnnData` object is not uploaded on the Github repo, but will be available alongside source code submission.  

Within source code submission, there will be 2 additional .h5ad objects in the folder:  
(1) `stephenson_2021_subsampled.h5ad`: this is equivalent to what users get from `pertpy` (downsampled data, no analysis)
(2) `adata_processed.h5ad`: this is processed `adata` object from running `01_preprocessing.ipynb`. This object is used in downstream ML models and contains representations, HVG expressions.