AutoGeneS
=========

AutoGeneS is a tool to automatically select informative genes from RNA-seq data [?]. Using this gene selection, it can be used to perform bulk deconvolution.

AutoGeneS requires no prior knowledge about marker genes and selects genes by simultaneously optimizing multiple criteria: minimizing the correlation and maximizing the distance between cell types [isn't this genes?]. It can be applied to reference profiles from various sources like single-cell experiments or sorted cell populations.

It has been designed to be compatible with `scanpy`_. To report issues or view the code, please refer to our `github`_ page.

.. _scanpy: https://github.com/theislab/scanpy
.. _github: https://github.com/theislab/AutoGeneS

.. toctree::
  :maxdepth: 1
  :hidden:

  background 
  getting_started 
  api
  references 

.. toctree::
  :maxdepth: 1
  :hidden:
  :caption: Tutorials
  
  tutorials/preprocessing
  tutorials/optimization
  tutorials/plotting_selection
  tutorials/deconvolution

.. toctree::
  :maxdepth: 1
  :hidden:
  :caption: Applications
  
  applications/bulk_deconvolution
