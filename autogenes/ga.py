import numpy as np
import random
import hashlib

from deap import creator, base, tools, algorithms
from cachetools import LRUCache

creator.create("FitnessGA", base.Fitness, weights=())
creator.create("IndividualGA", np.ndarray, fitness=creator.FitnessGA)
creator.IndividualGA.__eq__ = lambda self, other: np.array_equal(self,other)
creator.IndividualGA.__hash__ = lambda self: int(hashlib.sha1(self.view(np.uint8)).hexdigest(),16)

class GeneticAlgorithm:

  # Algorithm parameters
  PARAMETERS = {
      'population_size': 100,
      'offspring_size': 100,
      'crossover_pb': 0.7,
      'mutation_pb': 0.3,
      'mutate_flip_pb': 1E-3,
      'crossover_thres': 1000,
      'ind_standard_pb': 0.1
  }

  # Maximal number of individuals stored in cache
  MAX_CACHE_SIZE = 2E7

  def __init__(self, **kwargs):
    
    self.data = kwargs['data']
    self.ind_size = self.data.shape[1]
    self.min_nfeat = self.data.shape[0]

    self.objectives_func = kwargs['objectives_func']
    self.objectives_names = kwargs['objectives_names']

    self.ngen = kwargs['ngen']
    self.verbose = kwargs['verbose']

    creator.FitnessGA.weights = kwargs['weights']

    self._ncx = 0

    # Parameters
    self.params = GeneticAlgorithm.PARAMETERS.copy()
    for k in self.params:
      if k in kwargs:
        self.params[k] = kwargs[k]

    # Seed
    self.gen = np.random.RandomState(kwargs['seed'])
    # DEAP is using the python random module
    random.seed(kwargs['seed'])

    if not hasattr(creator, 'FitnessGA'):
      creator.create("FitnessGA", base.Fitness, weights=kwargs['weights'])
    else:
      creator.FitnessGA.weights = kwargs['weights']

    if not hasattr(creator, 'IndividualGA'):
      creator.create("IndividualGA", np.ndarray, fitness=creator.FitnessGA)
      creator.IndividualGA.__eq__ = lambda self, other: np.array_equal(self,other)
      creator.IndividualGA.__hash__ = lambda self: int(hashlib.sha1(self.view(np.uint8)).hexdigest(),16)

    self.tb = base.Toolbox()
    self.tb.register("evaluate", self.fitness)

    if kwargs['mode'] == 'standard':
      self.tb.register("individual", self.individual_standard, creator.IndividualGA)
      self.tb.register("mate", self.crossover_standard)
      self.tb.register("mutate", self.mutate_standard)
    elif kwargs['mode'] == 'fixed':
      self.tb.register("individual", self.individual_fixed, creator.IndividualGA)
      self.tb.register("mate", self.crossover_fixed)
      self.tb.register("mutate", self.mutate_fixed)
      self.nfeatures = kwargs['nfeatures']

    self.tb.register("select", tools.selNSGA2)

    # Initial population
    self.tb.register("population", tools.initRepeat, list, self.tb.individual)
    self.pop = self.tb.population(n=self.params['population_size'])

    # Pareto front
    self.hof = tools.ParetoFront(similar=np.array_equal)

    # Statistics
    self.stats = None
    self.count_gen = 0
    self.last_random_state = []

    self.stats = tools.Statistics(lambda ind: ind.fitness.values)
    
    # Prints out length of pareto and stores state for resume()
    def pareto_save_state(data):
      self.count_gen += 1
      self.last_ncx = self._ncx
      self.last_random_state = [self.gen.get_state(),random.getstate()]
      return str(len(self.hof.items))

    self.stats.register('pareto', pareto_save_state)

    for i,name in enumerate(self.objectives_names):
      def create_stats_func(i):
        def stats_func(data):
          data_i = list(map(lambda y: y[i], data))
          data_min = np.round(np.min(data_i),2)
          data_max = np.round(np.max(data_i),2)
          return f'{data_min} - {data_max}'
        return stats_func

      self.stats.register(name, create_stats_func(i))

    # Cache for fitness values
    self.tabu = LRUCache(GeneticAlgorithm.MAX_CACHE_SIZE)

  def run(self):

    try:
      self.pop, self.log = algorithms.eaMuPlusLambda(self.pop, self.tb,
                                   mu = self.params['population_size'], 
                                   lambda_ = self.params['offspring_size'],
                                   cxpb = self.params['crossover_pb'],
                                   mutpb = self.params['mutation_pb'],
                                   ngen = self.ngen - self.count_gen,
                                   stats = self.stats,
                                   verbose = self.verbose,
                                   halloffame = self.hof)
    except KeyboardInterrupt:
      # The first evaluation is the initial population
      self.count_gen -= 1
      if self.verbose:
        print()
        print("Stopped manually. Resume with ag.resume()")

    return self.hof

  def resume(self):

    print(f"Resume at generation {self.count_gen}")
    print()

    # Restore state
    self._ncx = self.last_ncx
    self.gen.set_state(self.last_random_state[0])
    random.setstate(self.last_random_state[1])

    self.run()

  def fitness(self, ind):

    selection = self.data[:,ind] 
    if ind not in self.tabu:
      obj = [f(selection) for f in self.objectives_func]
      self.tabu[ind] = obj
      return tuple(obj)
    else:
      return tuple(self.tabu[ind])

  def individual_standard(self, *args, **kwargs):
    """
    generate a random individual
    the number of features is binomially distributed, with at least self.min_nfeat features
    the features are then selected at random
    """

    ind = creator.IndividualGA(np.full(self.ind_size,False,dtype=bool))
    rate = self.params['ind_standard_pb']
    nfeat = self.min_nfeat + self.gen.binomial(self.ind_size-self.min_nfeat,rate)
    ind[:nfeat] = True
    self.gen.shuffle(ind)
    
    return ind
  
  def mutate_standard(self, ind):
    """
    mutate an individual by bit-flipping one or more randomly chosen
    elements
    """

    self.n_bitflip(ind)
    self.fill_up(ind)

    return ind,

  def crossover_standard(self, ind1, ind2):
    """
    mates two individuals
    either uses logical "and" and "or" or
    randomly swaps a block of rows or columns between two individuals
    """

    self._ncx += 1
    if self._ncx < self.params['crossover_thres']:
      self.cross_and_or(ind1, ind2)
    else:
      self.swap_block(ind1, ind2)

    self.fill_up(ind1)
    self.fill_up(ind2)

    return ind1, ind2

  def individual_fixed(self, *args, **kwargs):

    ind = creator.IndividualGA(np.full(self.ind_size,False,dtype=bool))
    ind[:self.nfeatures] = True
    self.gen.shuffle(ind)

    return ind

  def mutate_fixed(self, ind):

    # Calculate number of bits to flip
    # Later, n_flip 1s and n_flip 0s are flipped
    n_flip = self.gen.binomial(self.ind_size, self.params['mutate_flip_pb'])
    n_flip = int((n_flip-n_flip%2)/2)
    n_flip = min(self.nfeatures,self.ind_size-self.nfeatures,n_flip)

    # Choose which bits to flip
    id_1 = np.nonzero( ind )[0]
    id_0 = np.nonzero( np.invert(ind) )[0]

    ch_1 = self.gen.choice(id_1, n_flip, replace=False)
    ch_0 = self.gen.choice(id_0, n_flip, replace=False)

    # Create bitmask and apply it
    bm = np.array([False]*self.ind_size)
    bm[ch_1] = True
    bm[ch_0] = True

    np.logical_xor(ind, bm, out=ind)

    return ind,
  
  def crossover_fixed(self, ind1, ind2):

    union = np.logical_or(ind1, ind2)
    id_1 = np.nonzero(union)[0]
    ind1.fill(False)
    ind2.fill(False)
    ind1[id_1[:self.nfeatures]] = True
    ind2[id_1[-self.nfeatures:]] = True

    return ind1, ind2

  #
  # Helper
  #

  def swap_block(self, ind1, ind2):

    cx1 = self.gen.randint(0, self.ind_size)
    cx2 = self.gen.randint(cx1 + 1, self.ind_size + 1)

    ind1[cx1:cx2], ind2[cx1:cx2] = (
        ind2[cx1:cx2].copy(), ind1[cx1:cx2].copy())

  def cross_and_or(self, ind1, ind2):

    tmp = ind1.copy()
    np.logical_and(ind1,ind2,out=ind1)
    np.logical_or(ind1,ind2,out=ind1)

  def n_bitflip(self, ind):

    n_flip = self.gen.binomial(self.ind_size, self.params['mutate_flip_pb'])

    # generate bitmask
    bm = np.array([False]*self.ind_size)
    bm[1:n_flip] = True
    self.gen.shuffle(bm)
    
    # apply bitmask
    np.logical_xor(ind, bm, out=ind)

  def fill_up(self, ind):
    """
    fill up individual until it has at least min_nfeat features
    """

    add = self.min_nfeat - sum(ind)
    if add > 0:
      id_0 = np.nonzero(np.invert(ind))[0]
      ind[self.gen.choice(id_0, add, replace=False)] = True
