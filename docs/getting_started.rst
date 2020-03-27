Getting Started
===============

Install AutoGeneS with::

  pip install --user autogenes

In this section, we showcase the typical way our tool is used. We recommend running the code in your own Jupyter Notebook.

Start by reading in some sample data (available from our repository)::

  import anndata
  import numpy as np
  import pandas as pd
  import autogenes as amport anndata

  data = pd.read_csv('../datasets/GSE75748_bulk_data.csv',index_col='index')
  data = data.T.iloc[:,:100]
  adata = anndata.AnnData(data)
  adata.obs['celltype'] = adata.obs.index

Before you can use AutoGeneS, it needs to be initialized with the data. At this step, preprocessing happens (see the documentation)::

  ag.init(adata)

In the next step, we run the optimizer::

  ag.optimize(ngen=5)

Here, we only run it for five generation. During this optimization process, a set of solutions is generated. Each of them is evaluated based on objectives that can be passed to `optimize`. In our examples, it uses the default objectives `correlation` and `distance`.

Now, we want to plot the solution set::
  
  ag.plot(weights=(-1,1))

This plots the objective values of all solutions. The arguments are used to select a solution, which is marked in the plot. In our case, we choose the solution  using the weights `(-1,1)` on the objective values.

To choose a solution, run::
  
  ag.select(close_to=(1,75))

The following criterion is used to select a solution: The first number, `1` refers to the second objective, which is `distance` in our case. So, the above will choose the solution whose distance value is *closest* to 75.

Now that we have chosen a gene selection, we can deconvolve bulk data::

  coef = ag.deconvolve(adata.X[1], model='linear')
  print(coef[0])

This deconvolves a bulk sample that only contains the second cell type -- and it uses a linear model. The result, as predicted, is ``[0,1,0, ..., 0]``.

It is possible to compress all of the above steps into a single line::

  ag.pipeline(adata,adata.X[1],ngen=5,close_to=(1,75),model='linear')

This should produce the same result!

For more information on each step, please refer to the tutorial sections. For more extensive examples, see the examples section.
