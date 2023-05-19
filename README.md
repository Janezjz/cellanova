# CellANOVA: Cell State Space Analysis of Variance

- [Overview](#overview)
- [Installation](#installation)
- [Tutorials](#tutorials)
- [License](#license)

<img src="https://github.com/Janezjz/cellanova/blob/main/figures/model.jpg" width="700" height="450">

## Overview
The integration of cells across samples to remove unwanted batch variation plays a critical role in single cell analyses. When the samples are expected to be biologically distinct, it is often unclear how aggressively the cells should be aligned across samples to achieve uniformity.  CellANOVA is a Python package for batch integration with signal recovery in single cell data.  It builds on existing single cell data integration methods, and uses a pool of control samples to quantify the batch effect and separate meaningful biological variation from unwanted batch variation.  When used with an existing integration method, CellAnova allows the recovery of biological signals that are lost during integration.  

There are two ways to use CellANOVA.  You can start from scratch, or you can start with an existing integration.  If you start from scratch, CellANOVA will compute an initial integration with [Harmony](https://portals.broadinstitute.org/harmony/).  However, we have also achieved good results when starting with integration computed using [Seurat](https://satijalab.org/seurat/articles/integration_rpca.html).  The method is agnostic to the initial integration algorithm, and if you prefer to start by performing your own integration, you can choose any algorithm  hat gives a reasonable initial alignment of your data.

You will also need to select a pool of control samples.   These ``control'' samples will be used to estimate a latent linear space that captures cell- and gene-specific unwanted variations.  The basis vectors of this linear space will then be used to quantify the batch effect in each cell across all samples, and this will allow you to recover any cell- and gene-specific biological signal that is separable from batch effect that may have been erased during the initial integration.  By using only samples in the control pool in the estimation of the batch variation space, CellANOVA preserves any biological differences in the non-control samples that lie outside this space.  Importantly, CellANOVA produces a batch corrected gene expression matrix which can be used for gene- and pathway-level downstream analyses.
The selection of control pool samples depends on your experimental design.  Please see our manuscript for examples on how controls could be selected in case-control, longitudinal, or other single cell study designs.  

CellANOVA has been tested on data sets consisting of hundreds of thousands of cells.  Here are the runtimes recorded on a XX machine:
    100,000 cells:
    1,000,000 cells:

For more model details, validation results and real dataset analysis, please check out our paper (to be added). If you use our method, please use the following citation (to be added):

## Installation

Our Python package has been tested on python=3.7, 3.8, 3.9. It depends on numpy>=1.20.3, scipy>=1.7.1, pandas>=1.3.2, scikit-learn>=1.0.2, anndata>=0.7.6, scanpy>=1.8.1, harmonypy>=0.0.6. 

If you use conda environment, you can use the following command for an easy setup. It will build a seperate cellanova environment, and have all dependencies installed.

```bash
conda env create -f environment.yml
```

## Tutorials
### Quick Start

The following is a quick example showing CellANOVA integration pipeline.

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


### Example Workflows

For more detailed illustrations, please refer to the following jupyter notebooks:
 * [CellANOVA integration workfolow](https://github.com/Janezjz/cellanova/blob/main/tutorials/cellanova_integration.ipynb)
 * [Evaluation of batch removal performance](https://github.com/Janezjz/cellanova/blob/main/tutorials/eval_batch_removal.ipynb)
 * [Evaluation of distortion](https://github.com/Janezjz/cellanova/blob/main/tutorials/eval_distortion.ipynb)
 * [Evaluation of signal preservation](https://github.com/Janezjz/cellanova/blob/main/tutorials/eval_signal_preservation.ipynb)



## License

This project is licensed under the terms of GNU GENERAL PUBLIC LICENSE Version 3.
