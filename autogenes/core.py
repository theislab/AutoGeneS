import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from .ga import GeneticAlgorithm

from . import objectives as ga_objectives

import deap
import warnings

class AutoGeneS:

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

  def __init__(self, data):

    self.data = data

    if len(self.data.shape) != 2:
      raise ValueError("data is expected to have two dimensions")

    if self.data.shape[0] < 2:
      raise ValueError("At least two rows (cell types) expected")

    if self.data.shape[1] < self.data.shape[0]:
      raise ValueError("Number of columns (genes) must be >= number of rows (cell types)")

    if not np.isfinite(self.data).all():
      raise ValueError("Some entries in data are not scalars")

    self.__has_run = False
    self.selection = None
    self.selection_index = None

  def run(self, ngen=2, mode='standard', nfeatures=None, weights=None, objectives=None, seed=0, verbose=True, **kwargs):

    # Check modes

    if mode == 'standard':
      if nfeatures is not None:
        raise ValueError("nfeatures doesn't apply to standard mode (did you mean mode='fixed'?)")

    elif mode == 'fixed':
      if nfeatures is None:
        raise ValueError("You need to supply nfeatures")

      if nfeatures > self.data.shape[1]:
        raise ValueError("nfeatures must be <= the number of columns (genes)")
      
      if nfeatures < self.data.shape[0]:
        raise ValueError("nfeatures must be >= the number of rows (cell types)")
    else:
      raise ValueError("Invalid mode")

    # Check weights and objectives
    
    if weights is None:
      if objectives is None:
        weights = (-1.0,1.0)
        objectives = ('correlation','distance')
      else:
        raise Exception("Need weights for objectives")
    else:
      if objectives is not None:
        if len(weights) != len(objectives):
          raise ValueError("Number of weights does not match number of objectives")
        weights_l = []
        objectives_l = []
        for i,w in enumerate(weights):
          if w == 0:
            warnings.warn(f"Ignoring objective '{str(objectives[i])}'")
          else:
            weights_l.append(w)
            objectives_l.append(objectives[i])
        weights=tuple(weights_l)
        objectives=tuple(objectives_l)
      else:
        raise Exception("Need objectives for weights")
    
    # Store objectives

    self.objectives_func = []
    self.objectives_names = []

    for f in objectives:
      if callable(f):
        self.objectives_func.append(f)
        self.objectives_names.append(f.__name__)
      elif isinstance(f,str):
        if not hasattr(ga_objectives,f):
          raise ValueError(f"No such objective: {f}")
        else:
          self.objectives_names.append(f)
          self.objectives_func.append(getattr(ga_objectives,f))
      else:
        raise ValueError("Invalid objective")

    self.objectives_num = len(self.objectives_func)
    self.weights = weights

    self.ga = GeneticAlgorithm(
        data=self.data, 
        ngen=ngen,
        mode=mode,
        weights=weights, 
        objectives_names=self.objectives_names, 
        objectives_func=self.objectives_func, 
        seed=seed, 
        verbose=verbose,
        nfeatures=nfeatures,
        **kwargs
      )
    self.hof = self.ga.run()

    self.__has_run = True

  def resume(self):
    self.ga.resume()

  @property
  def pareto(self):
    self.__assert_run()
    return self.hof.items

  @property
  def fitness_matrix(self):
    self.__assert_run()

    all = []
    for i in range(self.objectives_num):
      vals = np.array(list(map(lambda x: x.fitness.values[i], self.hof.items)))
      all.append(vals)
    return np.array(all).T

  #
  # Plot results
  #

  def plot(self,objectives=(0,1), **kwargs):

    self.__assert_run()

    if self.objectives_num == 1:
      raise Exception("Cannot plot for a single objective")

    obj = objectives

    if len(obj) != 2:
      raise ValueError("Must supply two objectives per plot") 

    if not all(map(lambda x: x in range(self.objectives_num), obj)):
      raise ValueError(f"Invalid objectives, must be 0 <= x <= {self.objectives_num-1}")

    if not kwargs:
      return self.plot(weights=self.weights)

    i,desc = self.__from_pareto(**kwargs)

    if desc == 'index': legend = f'By index'
    if desc == 'weights': legend = f"Using weights {kwargs['weights']}"
    if desc == 'close_to': legend = f"Close to {kwargs['close_to'][1]}"

    if 'size' in kwargs:
      if kwargs['size'] not in ['small','large']:
        raise ValueError("Invalid size")
      size = kwargs['size']
    else:
      if len(self.pareto) < AutoGeneS.PLOT_THRESHOLD:
        size = 'small' 
      else:
        size = 'large'

    df = pd.DataFrame(self.fitness_matrix).sort_values(by=obj[0])

    df_0 = df[obj[0]]
    df_1 = df[obj[1]]

    params = AutoGeneS.PLOT_PARAMS[size]

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

  def select(self, **kwargs):
    self.__assert_run()

    if not kwargs:
      return self.select(weights=self.weights)

    i,desc = self.__from_pareto(**kwargs)
    self.selection = self.hof[i]
    self.selection_index = i

    return self.selection

  #
  # Helper
  #

  def __from_pareto(self,**kwargs):

    if sum([ x in kwargs for x in ["weights","index","close_to"]]) != 1:
      raise Exception("You need to provide exactly one criterion.")

    if 'weights' in kwargs:
      weights = kwargs['weights']
      i_max = self.__index_by_weights(weights) 
      return i_max,'weights'

    if 'index' in kwargs:
      index = kwargs['index']
      if isinstance(index,int):
        if index not in range(len(self.pareto)):
          raise ValueError("Invalid index")
        return index,'index'
      else:
        obj,i = index
        fit = pd.DataFrame(data=self.fitness_matrix).sort_values(by=obj)
        return fit.index.values[i],'index'
    
    if 'close_to' in kwargs:
      obj,num = kwargs['close_to']
      fit = self.fitness_matrix[:,obj]
      i = np.argmin(np.abs(fit-num))
      return i,'close_to'

  def __index_by_weights(self,weights):
    self.__assert_run()
    if len(weights) != self.objectives_num:
      raise ValueError(f"Number of weights does not match number of objectives")

    fitness = self.fitness_matrix
    for i in range(self.objectives_num):
      max = np.max(fitness[:,i])
      if max:
        fitness[:,i] *= 1/max

    wfitness = fitness.dot(np.array(weights))
    return np.argmax(wfitness)

  def __assert_run(self):
    if not self.__has_run:
      raise Exception("AutoGeneS did not run yet")

  def __setstate__(self,dict):
    deap.creator.FitnessGA.weights = dict['weights']
    self.__dict__.update(dict)
