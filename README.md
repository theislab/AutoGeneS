# AutoGeneS

AutoGeneS **automatically** extracts informative genes and reveals the cellular heterogeneity of bulk RNA samples. AutoGeneS requires no prior knowledge about marker genes and selects genes by **simultaneously optimizing multiple criteria**: minimizing the correlation and maximizing the distance between cell types. It can be applied to reference profiles from various sources like single-cell experiments or sorted cell populations.

For a multi-objective optimization problem, there usually exists no single solution that simultaneously optimizes all objectives. In this case, the objective functions are said to be conflicting, and there exists a (possibly infinite) number of **Pareto-optimal solutions**. Pareto-(semi)optimal solutions are a set of all solutions that are not dominated by any other explored solution. Pareto-optimal solutions offer a set of equally good solutions from which to select, depending on the dataset

## Dependencies

* python >=3.6
* dill
* deap
* cachetools

## Installation
git clone https://github.com/theislab/AutoGeneS<br/>
pip install --user dist/autogenes-0.9.1-py3-none-any.whl<br/>

## Testing
```python
import numpy as np
from autogenes import AutoGenes
ag = AutoGenes(np.identity(2))
```

## Usage
The normalized reference profiles filtered for highly variable genes are given as input to AutoGeneS as numpy array or panda dataframe with the format: cell_types x genes.
We recommend to perform the optimization on 4,000-5,000 highly variable genes.

```python
from autogenes import AutoGenes
ag = AutoGenes(centroids_hv)
ag.run(ngen=5000,seed=0,nfeatures=400,mode='fixed') #ngen is the number of optimization runs and nfeatures is the number of marker genes we are interested in
```

The pareto front solutions are then accessible as follows:<br/>
```python
pareto = ag.pareto
``` 

We then pick one solution and filter its corresponding marker genes:<br/>
```python
centroids_pareto = centroids_hv.T[pareto[len(pareto)-1]] #select the solution with min correlation
``` 

For deconvolution, bulk RNA samples are first normalized and filtered for the same marker genes. The deconvolution is then performed as follows:
```python
regr_NuSVR = NuSVR(nu=0.5,C=0.5,kernel='linear') 
regr_NuSVR.fit(centroids_pareto, bulk_pareto)
``` 

Later, ```regr_NuSVR.coef_[0]``` returns the proportions which should be normalized to sum to one. 

## Example Notebooks

Simple examples for how to call AutoGeneS can be find in "jupyter"

Bulk deconvolution of synthetic bulks using single-cell reference profiles can be find in "deconv_example"
