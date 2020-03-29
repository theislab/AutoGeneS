Getting Started
===============

Install AutoGeneS with::

  pip install --user autogenes

In the following, we show how to use AutoGeneS with an example. 

Import the libraries and read the reference data and bulk samples::

  import anndata
  import numpy as np
  import pandas as pd
  import autogenes as ag
  import anndata

  bulk_data = pd.read_csv('bulk_data.csv').transpose()
  adata = sc.read(address_to_your_sc_data, cache=True).transpose()

Before you can use AutoGeneS, it needs to be initialized with the reference data (see the API)::

  ag.init(adata,use_highly_variable=True,celltype_key='cellType')

If the data is given as anndata, ag.init automatically measures the centroids of cell types by means of averaging the gene expression of their cells.

In the next step, we run the optimizer::

  ag.optimize(ngen=5000,nfeatures=400,seed=0,mode='fixed')

Here, we run the optimizer for 5K generations asking for 400 genes. During this optimization process, a set of solutions is generated. Each solution is a set of 400 genes and is evaluated based on objectives that can be passed to `optimize`. In our examples, it uses the default objectives `correlation` and `distance`.

All the non-dominated solutions can be visualized using the plot function::
  
  ag.plot(weights=(-1,0))

This plots the objective values of all non-dominated solutions. The arguments are used to select a solution, which is marked in the plot. In our case, we choose the solution  using the weights `(-1,0)` on the objective values that will return the solution with minimum correlation.

To choose another solution, run::
  
  ag.select(close_to=(1,75))

The following criterion is used to select a solution: The first number, `1` refers to the second objective, which is `distance` in our case. So, the above will choose the solution whose distance value is *closest* to 75. There are other ways of selecting a solution like specifying the index of a solution using `ag.select(index=0)`.

Now that we have chosen a solution with a set of genes, we can deconvolute bulk samples::

  coef = ag.deconvolve(bulk_data, model='nnls')
  print(coef)

We recommend to normalize the coeffiecients after the analysis.

It is possible to compress all of the above steps into a single line::

  ag.pipeline(adata,bulk_data,ngen=5000,nfeatures=400,seed=0,mode='fixed',close_to=(1,75),model='nnls')

This should produce the same result!

For more information on each step, please refer to the API. For more extensive examples, see the examples section.
