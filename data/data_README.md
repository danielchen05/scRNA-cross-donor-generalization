This project uses datasets accessed through the `pertpy` package through the commands:  
```python
import pertpy as pt
adata = pt.data.stephenson_2021_subsampled()
# adata.write("../data/stephenson_2021_subsampled.h5ad")
```
Data could be read in in two ways. We will provide code to read in data in both ways in each notebook/script, with default being option 1 (using `pertpy`).   
(1) Using the commands above.   
(2) Raw data is also saved as a `AnnData` object to avoid dependency of `pertpy` package (it works weird on Windows laptops). The `AnnData` object is not uploaded on the Github repo, but will be available alongside source code submission.  