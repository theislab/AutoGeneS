from .core import AutoGeneS
from typing import Optional, Tuple

import pandas as pd
import anndata
import numpy as np
import warnings
import dill as pickle

from sklearn.svm import NuSVR
from sklearn import linear_model
from scipy.optimize import nnls

class Interface:

  def __init__(self):

    self.pre_selection = None
    self._selection = None

    self._adata = None
    self.data = None
    self.data_genes = None

    self.main = None

  def init(self, 
      data, 
      celltype_key = 'celltype', 
      genes_key = None,
      use_highly_variable = False, 
      **kwargs
    ):
    """
    init(data, celltype_key = 'celltype', genes_key = None, use_highly_variable = False)
    
    Preprocesses input data

    If an AnnData object is passed, it is assumed that it contains single-cell data. The means are calculated using 'celltype_key'. In addition, a pre-selection of genes can be specified with 'genes_key' or 'use_highly_variable'. Then, only these genes will be considered in the optimization.

    If a DataFrame or numpy array is passed, it is assumed that they already contain the means.

    Parameters
    ----------
    data : `anndata.AnnData`, `np.ndarray`, `pd.DataFrame`
      Input data
    celltype_key : `str`, optional (default: `celltype`)
      Name of the obs column that specifies the cell type of a cell
      For AnnData only
    genes_key : `str`, optional (default: `None`)
      Name of the var column with boolean values to pre-select genes 
    use_highly_variable : `bool`, optional (default: `False`)
      Equivalent to genes_key='highly_variable'

    Returns
    -------
    None
    """

    self.__init__()
    
    # Process different input formats
    if isinstance(data, anndata.AnnData):

      if use_highly_variable: genes_key = 'highly_variable'

      if celltype_key not in data.obs:
        raise ValueError(f"AnnData has no obs column '{celltype_key}'")

      self._adata = self.__compute_means(data,celltype_key)
      self.data_genes = data.var_names.values

      if genes_key:
        self.pre_selection = data.var[genes_key].values
      else:
        self.pre_selection = np.full((data.X.shape[1],),True)

      self.data = self._adata.X[:,self.pre_selection]

      self.main = AutoGeneS(self.data)
      return self._adata

    elif isinstance(data,pd.DataFrame):

      self.data = data.values
      self.data_genes = data.columns.values
      self.main = AutoGeneS(self.data)
      self.pre_selection = np.full((data.shape[1],),True)

    elif isinstance(data, np.ndarray):
      
      self.data = data
      self.main = AutoGeneS(self.data)
      self.pre_selection = np.full((data.shape[1],),True)

    else:
      raise TypeError("data must be AnnData, DataFrame or ndarray")

  def optimize(
      self, 
      ngen = 2,
      mode = 'standard',
      nfeatures = None,
      weights = None,
      objectives = None,
      seed = 0,
      verbose = True,
      **kwargs
  ):
    """
    optimize(ngen = 2, mode = 'standard', nfeatures = None, weights = None, objectives = None, seed = 0, verbose = True, **kwargs)

    Runs multi-objective optimizer

    This method runs an evolutionary algorithm to find gene selections that optimize certain objectives. It can run for a different number of generations and in different modes. For more information on genetic algorithms and their parameters, refer to the `DEAP documention <https://deap.readthedocs.io/en/master/index.html>`_.

    Parameters
    ----------
    ngen : `int`, optional (default: `2`)
      Number of generations. The higher, the longer it takes
    mode : `standard`, `fixed`, optional (default: `standard`)
      In standard mode, the number of genes of a selection is allowed to vary arbitrarily. In fixed mode, the number of selected genes is fixed (using `nfeatures`)
    nfeatures : `int`, optional (default: `int`)
      Number of genes to be selected in fixed mode
    weights : `(int, ...)`, optional (default: `(-1,1)`)
      Weights applied to the objectives. For the optimization, only the sign is relevant: `1` means to maximize the respective objective, `-1` to minimize it and `0` means to ignore it. The weight supplied here will be the default weight for selection. There must be as many weights as there are objectives
    objectives : `([str,function], ...)`, optional (default: `('correlation','distance')`)
      The objectives to maximize or minimize. Must have the same length as weights. The default objectives (correlation, distance) can be referred to using strings. For custom objectives, a function has to be passed. For further details, refer to the respective tutorial.
    seed : `int`, optional (default: `0`)
      Seed for random number generators
    verbose : `bool`, optional (default: `True`)
      If True, output a progress summary of the optimization (the current generation, size of the pareto front, min and max values of all objectives)
    population_size : `int`, optional (default: `100`)
      Size of every generation (mu parameter)
    offspring_size : `int`, optional (default: `50`)
      Number of individuals created in every generation (lambda parameter)
    crossover_pb : `float`, optional (default: `0.7`)
      Crossover probability
    mutation_pb : `float`, optional (default: `0.3`)
      Mutation probability
    mutate_flip_pb : `float`, optional (default: `1E-3`)
      Mutation flipping probability (fixed mode)
    crossover_thres : `int`, optional (default: `1000`)
      Crossover threshold (standard mode) 
    ind_standard_pb : `float`, optional (default: `0.1`)
      Probability used to generate initial population in standard mode

    Returns
    -------
    None
    """

    if self.main is None:
      raise Exception("Not initialized")

    self.main.run(
        ngen=ngen,
        mode=mode,
        nfeatures=nfeatures,
        weights=weights,
        objectives=objectives,
        seed=seed,
        verbose=verbose,
        **kwargs
      )

  def plot(self, **kwargs):
    """
    plot(objectives = (0,1), weights = None, index = None, close_to = None)
    
    Plots objective values of solutions

    Can only be run after `optimize`. Every parameter corresponds to one selection method. Only one can be chosen at a time. If you don't specify an selection method, the weights passed to `optimize` will be used.

    Parameters
    ----------
    objectives : `(int,int)`, optional (default: `(0,1)`)
      The objectives to be plotted. Contains indices of objectives. The first index refers to the objective that is plotted on the x-axis. For example, `(2,1)` will plot the third objective on the x-axis and the second on the y-axis.
    weights : `(int, ...)`, optional
      Weights with which to weight the objective values. For example, `(-1,2)` will minimize the first objective and maximize the the second (with higher weight).
    index : `int`, `(int,int)`, optional
      If one int is passed, return `pareto[index]`
      If two ints are passed, the first is an objective (`0` for the first). The second is the nth element if the solutions have been sorted by the objective in ascending order. For example, `(0,1)` will return the solution that has the second-lowest value in the first objective. `(1,-1)` will return the solution with the highest value in the second objective.
    close_to : `(int,int)`, optional
      Select the solution whose objective value is closest to a certain value. Assumes `(objective,value)`. For example, `(0,100)` will select the solution whose value for the first objective is closest to 100.
    """

    if self.main is None:
      raise Exception("Not initialized")

    self.main.plot(**kwargs)

  def select(self, copy=False, key_added='autogenes', **kwargs):
    """
    select(weights = None, close_to = None, index = None, copy=False, key_added='autogenes')

    Selects a solution

    Specify a criterion to choose a solution from the solution set. Supports adding the solution to the annotation of an adata object. Can only be run after `optimize`

    Parameters
    ----------
    weights : `(int, ...)`, optional
      Weights with which to weight the objective values. For example, `(-1,2)` will minimize the first objective and maximize the the second (with more weight).
    index : `int`, `(int,int)`, optional
      If one int is passed, return `pareto[index]`
      If two ints are passed, the first is an objective (`0` for the first). The second is the nth element if the solutions have been sorted by the objective in ascending order. For example, `(0,1)` will return the solution that has the second-lowest value in the first objective. `(1,-1)` will return the solution with the highest value in the second objective.
    close_to : `(int,int)`, optional
      Select the solution whose objective value is close to a certain value. Assumes `(objective,value)`. For example, `(0,100)` will select the solution whose value for the first objective is closest to 100.
    copy : `bool`, optional (default: `False`)
      If true, a new adata object will be created with the selected solution in the var column specified by `key_added`
    key_added : `str`, optional (default: `autogenes`)
      The name of the var column to which to add the chosen gene selection
    """
    if self.main is None:
      raise Exception("Not initialized")

    s = self.main.select(**kwargs)

    self._selection = self.__process_selection(s)
    
    if self._adata:
      if copy:
        r = self._adata.copy()
        r.var[key_added] = self._selection
        return r
      else:
        self._adata.var[key_added] = self._selection

    return self._selection

  def deconvolve(self, bulk, key=None, model='nusvr', **kwargs):
    """
    deconvolve(bulk,key = None, model='nusvr')

    Performs bulk deconvolution

    Deconvolves bulk data using a gene selection. The selection can be specified through a key or the current selection is used.

    If the optimizer has been run, but nothing has been selected yet, an automatic selection occurs (equivalent to ``ag.select()``)

    Parameters
    ----------
    bulk : `np.ndarray`, `pd.Series`, `pd.DataFrame`, `AnnData`
      If multi-dimensional, then each row corresponds to a sample. If it has gene annotations (e.g. var_names for AnnData or df.columns for DataFrame), the method will respond intelligently (reorder if necessary, use only those genes from the selection that are available in the bulk data)
    key : `str`, optional (default: `None`)
      Name of the var column that specifies a gene selection. If None, then the current selection is used (or is automatically chosen)
    model : `nusvr`, `nnls`, `linear`, optional (default: `nusvr`)
      Choose a regression model. Available options: NuSVR, non-negative least squares and linear model.

    Returns
    -------
    An array of the form `[[float, ...],...]` containing the model coefficients for each target (bulk sample)
    """

    if self._selection is None:
      self.select(**kwargs)

    selection = self._adata.var[key] if key else self._selection

    bulk_data, bulk_genes = self.__unpack_bulk(bulk)

    X,y = self.__model_input(bulk_data, bulk_genes, selection)

    if model == "nusvr":
      if y.shape[1] == 1:
        y = np.ravel(y)
        model = NuSVR(nu=0.5,C=0.5,kernel='linear')
        model.fit(X, y)
        self.model = model
        return model.coef_
      else:
        res = np.zeros((y.shape[1],X.shape[1]))
        for i in range(y.shape[1]):
          model = NuSVR(nu=0.5,C=0.5,kernel='linear')
          model.fit(X, y[:,i])
          self.model = model
          res[i] = model.coef_
        return res

    if model == "nnls":
      if y.ndim == 1:
        x,err = nnls(X,y)
        return x
      else:
        res = np.zeros((y.shape[1],X.shape[1]))
        for i in range(y.shape[1]):
          x,err = nnls(X,y[:,i])
          res[i] = x
        return res

    if model == "linear":
      model = linear_model.LinearRegression(copy_X=True, fit_intercept=False)
      model.fit(X, y)
      self.model = model
      return model.coef_

    raise ValueError("Model is not supported")

  def pipeline(self, data, bulk, **kwargs):
    """
    pipeline(data,bulk, **kwargs)

    Runs the optimizer, selection and deconvolution using one method
    """

    self.init(data,**kwargs)
    self.run(**kwargs)
    return self.deconvolve(bulk, **kwargs)

  def resume(self):
    """Resumes an optimization process that has been interrupted"""

    if self.main is None:
      raise Exception("Not initialized")

    self.main.resume()

  def save(self,filename):
    """Saves current state to a file
    
    Parameters
    ----------
    filename : `str`
      Name of the file
    """
    pickle.dump(self, open(filename, 'wb'))

  def load(self,filename):
    """Loads a state from a file
    
    Parameters
    ----------
    filename : `str`
      Name of the file
    """
    tmp = pickle.load(open(filename, 'rb'))
    self.__dict__.update(tmp.__dict__)

  def adata(self):
    """Returns AnnData object
    
    Returns
    -------
    The AnnData object that the optimizer operates on (if no AnnData was passed to `ag.init`, `None`)
    """

    return self._adata

  def fitness_matrix(self):
    """Returns fitness matrix

    Returns
    -------
    A `pd.DataFrame` that contains the objective values of all solutions. The nth row corresponds to the nth solution (``ag.pareto()[n]``)
    """
    return self.main.fitness_matrix

  def pareto(self):
    """Returns the entire solution set
    
    Returns
    -------
    The solution set in the form `[[bool],...]`. Every member corresponds to a gene selection 
    """

    if self.main is None:
      raise Exception("Not initialized")

    return list(map(self.__process_selection, self.main.pareto))

  def selection(self):
    """Returns the current selection

    Returns
    -------
    The current selection as a boolean array
    """
    if self._selection is None:
      raise Exception("Nothing selected")
    return self._selection

  #
  # Helper
  #

  def __process_selection(self,s):

    r = self.pre_selection.copy()
    i = 0
    for k,val in enumerate(self.pre_selection):
      if val:
        r[k] = s[i]
        i += 1

    return r

  def __compute_means(self,adata,celltype_key):
    """
    returns a new, shallow (!) AnnData. It contains the mean gene expressions per cell type. The row names are the cell types. The column names are the genes of the original adata. 
    """

    if celltype_key not in adata.obs:
      raise ValueError("Key not found")
    sc_means = pd.DataFrame(data=adata.X, columns=adata.var_names)
    sc_means['cell_types'] = pd.Series(data=adata.obs[celltype_key].values,index=sc_means.index)
    sc_means = sc_means.groupby('cell_types').mean()

    if len(sc_means.index) == 1:
      raise ValueError("More than 1 cell types expected")

    result = anndata.AnnData(sc_means)
    result.var = adata.var.copy()
    result.var_names = adata.var_names

    return result

  def __model_input(self,bulk_data,bulk_genes,selection):

    data_genes = self.data_genes

    # Case: gene labels for both bulk and data are available
    if bulk_genes is not None and data_genes is not None:

      common_genes = np.isin(data_genes,bulk_genes)
      intersect_genes = np.logical_and(common_genes,selection)
      n_intersect_genes = sum(intersect_genes)

      if n_intersect_genes == 0:
        raise ValueError("None of the selected genes appear in the bulk data")

      if n_intersect_genes < sum(selection):
        warnings.warn("Some of the selected genes don't appear in the bulk data and will be ignored")

      if self._adata:
        X = self._adata.X.T[intersect_genes]
      else:
        X = self.data.T[intersect_genes]

      # Note: Genes in bulk may be in different order and of different size!
      # Cannot simply apply bitmask!

      y = np.zeros((bulk_data.shape[0],n_intersect_genes))
      gene_names = data_genes[intersect_genes]
      for i,gene in enumerate(gene_names):
        bulk_gene_index = np.argwhere(bulk_genes == gene)[0][0]
        y[:,i] = bulk_data[:, bulk_gene_index]

      y = y.T

    # Case: no gene labels available (for at least one)
    else:

      bulk_dim = bulk_data.shape[1]
      if bulk_dim != len(selection): #self.data.shape[1]
        raise ValueError("Bulk data has wrong shape")

      if self._adata:
        X = self._adata.X.T[selection]
      else:
        X = self.data.T[selection]
      y = bulk_data.T[selection]

    return X,y

  def __unpack_bulk(self,bulk):
    """
    returns tuple of
    2-dim ndarray bulk_data
    1-dim bulk_genes (or None)
    """

    bulk_data, bulk_genes = None, None

    if isinstance(bulk,np.ndarray):
      if bulk.ndim == 1:
        bulk = bulk.reshape(1,len(bulk))
      bulk_data = bulk

    if isinstance(bulk, anndata.AnnData):
      bulk_genes = bulk.var.index.values
      bulk_data = bulk.X

    if isinstance(bulk, pd.Series):
      bulk_genes = bulk.index.values
      bulk_data = bulk.values.reshape(1,len(bulk))

    if isinstance(bulk, pd.DataFrame):
      bulk_genes = bulk.columns.values
      bulk_data = bulk.values

    if bulk_data is None:
      raise ValueError("Invalid data type for bulk")

    return bulk_data,bulk_genes
