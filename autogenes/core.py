import pandas as pd
import numpy as np
#import anndata
import dill as pickle
import matplotlib.pyplot as plt

from .ga import GeneticAlgorithm

from . import objectives as ga_objectives

from sklearn.svm import NuSVR
from deap import creator

class AutoGenes:

  def __init__(self, data):

    self.adata = None
    self.dataframe = None

    #if isinstance(data, anndata.AnnData):
      #self.adata = data
      #self.data = data.X
    if isinstance(data,pd.DataFrame):
      self.dataframe = data
      self.data = data.values
    elif isinstance(data, np.ndarray):
      self.data = data
    else:
      raise TypeError("data must be AnnData, DataFrame or ndarray")

    if len(self.data.shape) != 2:
      raise ValueError("data is expected to have two dimensions")

    if self.data.shape[0] < 2:
      raise ValueError("At least two rows (cell types) expected")

    if self.data.shape[1] < self.data.shape[0]:
      raise ValueError("Number of columns (genes) must be >= number of rows (cell types)")

    if not np.isfinite(self.data).all():
      raise ValueError("Some entries in data are not scalars")

    self.__has_run = False
    self.__selection = None

  def run(self, nMU=100, ngen=2, mode='standard', nfeatures=None, weights=None, objectives=None, seed=0, verbose=True):

    if not isinstance(ngen, int) or ngen < 0:
      raise ValueError("Invalid ngen")

    if not isinstance(verbose, bool):
      raise ValueError("Invalid verbose")

    #
    # Modes
    #

    if mode == 'standard':
      if nfeatures is not None:
        raise ValueError("nfeatures doesn't apply to standard mode (did you mean mode='fixed'?)")

    elif mode == 'fixed':
      if nfeatures is None:
        raise ValueError("You need to supply nfeatures")

      if not isinstance(nfeatures, int):
        raise TypeError("nfeatures must be an integer")

      if nfeatures > self.data.shape[1]:
        raise ValueError("nfeatures must be <= the number of columns (genes)")
      
      if nfeatures < self.data.shape[0]:
        raise ValueError("nfeatures must be >= the number of rows (cell types)")
    else:
      raise ValueError("Invalid mode")

    #
    # Weights and objectives
    #
    
    if weights is None and objectives is not None:
      raise Exception("Need weights for objectives")

    if weights is not None and objectives is None:
      raise Exception("Need objectives for weights")
    
    if weights is None:
      weights = (1.0,-1.0)

    if objectives is None:
      objectives = ('distance','correlation')

    if (not isinstance(objectives, tuple)) or len(objectives) == 0:
      raise ValueError("objectives must be a nonempty tuple")
    
    self.objectives_func = []
    self.objectives_names = []

    for f in objectives:
      if callable(f):
        self.objectives_func.append(f)
        self.objectives_names.append(f.__name__)
      elif isinstance(f,str):
        if not hasattr(ga_objectives,f):
          raise ValueError(f"No such objective: {f}")
        self.objectives_names.append(f)
        self.objectives_func.append(getattr(ga_objectives,f))
      else:
        raise ValueError("Invalid objective")

    self.objectives_num = len(self.objectives_func)

    self.__check_weights(weights)
    self.weights = weights

    self.ga = GeneticAlgorithm(
        data=self.data, 
        nMU=nMU,
        ngen=ngen,
        mode=mode,
        weights=weights, 
        objectives_names=self.objectives_names, 
        objectives_func=self.objectives_func, 
        seed=seed, 
        verbose=verbose,
        nfeatures=nfeatures
      )
    self.hof = self.ga.run()

    self.__has_run = True

  def resume(self):
    self.ga.resume()

  #
  # Access results
  #

  @property
  def pareto(self):
    self.__assert_run()
    return self.hof.items

  def print(self):
    fitness = self.__fitness_matrix()

    for i in range(len(fitness)):
      val_str = '(' + ', '.join(fitness[i].astype(str)) + ')'
      print('%i: %s' %(i, val_str))

  def summary(self):
    fitness = self.__fitness_matrix()
    print(f'Pareto front contains {len(fitness)} individuals')
    for i,obj in enumerate(self.objectives_names):
      col = fitness[:,i]
      stats = np.round([np.min(col), np.mean(col), np.std(col), np.max(col)], 2)
      print(f"Objective '{obj}':")
      print(f"Min: {stats[0]}, Mean: {stats[1]}, Std: {stats[2]}, Max: {stats[3]}")

  #
  # Plot results
  #

  PLOT_PARAMS = {
        'small': {
            'figsize': (10,5),
            'all_ms': 8,
            'sel_ms': 10
        },
        'large': {
            'figsize': (15,10),
            'all_ms': 5,
            'sel_ms': 10
        }
    }

  PLOT_THRESHOLD = 50


  def plot(self,objectives=(0,1), **kwargs):

    self.__assert_run()

    if not isinstance(objectives, tuple):
      raise ValueError("objectives must be tuple")

    if len(objectives) != 2:
      raise ValueError("Must supply two objectives per plot") 

    if not all(map(lambda x: x in range(self.objectives_num), objectives)):
      raise ValueError(f"Invalid objectives, must be 0 <= x <= {self.objectives_num-1}")

    obj = objectives

    if 'index' in kwargs:
      i = kwargs['index']
      if i not in range(len(self.pareto)):
        raise ValueError("Invalid index")
      legend = f'Index {i}'
    else:
      weights = kwargs.get('weights', self.weights)
      self.__check_weights(weights)
      legend = f"Using weights {weights}"
      i = self.__index_by_weights(weights)

    if 'size' in kwargs:
      if kwargs['size'] not in ['small','large']:
        raise ValueError("Invalid size")
      size = kwargs['size']
    else:
      if len(self.pareto) < AutoGenes.PLOT_THRESHOLD:
        size = 'small' 
      else:
        size = 'large'

    df = pd.DataFrame(self.__fitness_matrix()).sort_values(by=obj[0])

    df_0 = df[obj[0]]
    df_1 = df[obj[1]]

    params = AutoGenes.PLOT_PARAMS[size]

    plt.figure(figsize=params['figsize'])

    line = plt.plot(df_0,df_1)

    plt_all, = plt.plot(df_0.drop(i),df_1.drop(i),'bo',ms=params['all_ms'])
    plt_sel, = plt.plot(df_0[i],df_1[i],'r^',ms=params['sel_ms'])

    plt.xlabel(self.objectives_names[obj[0]])
    plt.ylabel(self.objectives_names[obj[1]])

    plt.legend([plt_all, plt_sel], ["Option", legend],bbox_to_anchor=(1, 1), loc='upper left')

    plt.show()

  #
  # Select individual
  #

  def select(self, col='autogenes', **kwargs):
    self.__assert_run()

    if not kwargs:
      return self.select(weights=self.weights)

    if 'weights' in kwargs:
      weights = kwargs['weights']
      i_max = self.__index_by_weights(weights) 
      s = self.hof[i_max]
    elif 'index' in kwargs:
      i = kwargs['index']
      if i not in range(len(self.pareto)):
        raise ValueError("Invalid index")
      s = self.hof[i]
    else:
      raise ValueError("Invalid arguments")

    if self.adata:
      self.adata.var[col] = s
    self.__selection = s

    return s

  #
  # Bulk deconvolution
  #

  def bulk_deconvolution(self, bulk, model=None):
    """
    bulk is expected to be in the format that the model expects
    this is typically array_like and shape (n_genes, n_targets)
    It is not assumed that the selection has already been applied to bulk!
    """

    if self.__selection is None:
      raise Exception("Nothing selected yet")

    # Additional check if index was supplied both times
    if isinstance(bulk,pd.Series) and self.dataframe is not None:
      if any(bulk.index != self.dataframe.columns):
        raise ValueError("bulk is a DataFrame with different columns than dataframe")

    # Model
    if model is None:
      model = NuSVR(nu=0.5,C=0.5,kernel='linear')
    if not hasattr(model, 'fit'):
      raise ValueError("Invalid model")

    # Change to sklearn format
    X = self.data.T[self.__selection]
    y = bulk.T[self.__selection]

    model.fit(X, y)
    self.model = model
    return model.coef_

  def save(self, file):
    pickle.dump(self, open(file, 'wb'))

  def load(file):
    ag = pickle.load(open(file, 'rb'))
    creator.FitnessGA.weights = ag.weights
    return ag

  #
  # Helper
  #

  def __index_by_weights(self,weights):
    self.__assert_run()
    self.__check_weights(weights)

    fitness = self.__fitness_matrix()
    for i in range(self.objectives_num):
      max = np.max(fitness[:,i])
      if max:
        fitness[:,i] *= 1/max

    wfitness = fitness.dot(np.array(weights))
    return np.argmax(wfitness)

  def __fitness_matrix(self):
    self.__assert_run()

    all = []
    for i in range(self.objectives_num):
      vals = np.array(list(map(lambda x: x.fitness.values[i], self.hof.items)))
      all.append(vals)
    return np.array(all).T

  def __assert_run(self):
    if not self.__has_run:
      raise Exception("AutoGenes did not run yet")

  def __check_weights(self, weights):
    if not isinstance(weights, tuple):
      raise ValueError("weights must be a tuple")

    if len(weights) != self.objectives_num:
      raise ValueError(f"Wrong number of weights")

    if any([not (isinstance(w,(int,float)) or w == 0) for w in weights]):
      raise ValueError("Invalid weight value")
