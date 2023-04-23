# CellANOVA: Cell State Space Analysis of Variance

- [Introduction](#introduction)
- [Installation](#installation)
- [Tutorials](#tutorials)
- [License](#license)



## Introduction

CellANOVA is a Python package, for batch correction with signal recovery. It contructs a pool of control samples to estimate a latent linear space that captures cell- and gene-specific unwanted variations, which can then be used to remove batch effects from cells across all samples.  By using only samples in the control pool in the estimation of the batch variation space, CellANOVA preserves any biological differences in the non-control samples that lie outside this space.  Importantly, CellANOVA produces a batch corrected gene expression matrix which can be used for gene- and pathway-level downstream analyses, and is fast and scalable to data sets containing millions of cells. 

CellANOVA can be applied to multiple settings:
* Case-control design
* Longitudinal design
* Irregular block design

For more model details, validation results and real dataset analysis, please check out our paper (to add link). If you use our method, please use the following citation:




## Installation
### Dependencies

Our Python package has been tested on python=3.7, 3.8, 3.9. It depends on numpy>=1.20.3, scipy>=1.7.1, pandas>=1.3.2, scikit-learn>=1.0.2, anndata>=0.7.6, scanpy>=1.8.1, harmonypy>=0.0.6. 

If you use conda environment, you can use the following command for an easy setup. It will build a seperate cellanova environment, and have all dependencies installed.

```bash
conda env create -f environment.yml
```

### Download



## Tutorials
### Quick Start

The following is an quick example showing CellANOVA integration pipeline.

```python
## load required package
import anndata as ad
import scanpy as sc
from cellanova import *

## load and preprocess data
adata = sc.read_h5ad('raw_data.h5ad')
adata_prep = preprocess_data(adata, integrate_key='dataidx')

## construct control pool
control_dict = {
    'pool1': list(set(adata_prep[adata_prep.obs['condition']=='control',].obs['dataidx'])),
}

## model fitting
adata_prep= calc_ME(adata_prep, integrate_key='dataidx')
adata_prep = calc_BE(adata_prep, integrate_key, control_dict)
adata_prep = calc_TE(adata_prep, integrate_key)

## create an independent anndata object for cellanova-integrated data
integrated = ad.AnnData(adata_prep.layers['denoised'], dtype=np.float32)
integrated.obs = adata_prep.obs.copy()
integrated.var_names = adata_prep.var_names
```


### Example Notebooks

For more detailed examples, please refer to the following jupyter notebooks:
 * CellANOVA integration workfolow (link)
 * Evaluation of batch removal performance (link)
 * Evaluation of distortion (link)
 * Evaluation of signal preservation (link)



## License


