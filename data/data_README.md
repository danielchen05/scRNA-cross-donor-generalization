This project uses datasets accessed through the `pertpy` package through the commands:  
```python
import pertpy as pt
adata = pt.data.stephenson_2021_subsampled()
adata.write("../data/stephenson_2021_subsampled.h5ad")
```
Raw data is saved as a `AnnData` object to avoid future dependency of `pertpy` package. The `AnnData` object is not uploaded on the Github repo, but will be available alongside source code submission.