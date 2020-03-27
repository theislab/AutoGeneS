import pytest
import numpy as np
import pandas as pd

import os
import sys
sys.path.insert(0, '..')

from autogenes import AutoGeneS
from autogenes.ga import GeneticAlgorithm

@pytest.fixture
def ag():

  ag = AutoGeneS(np.identity(10)[:5])
  return ag

def test_ga_params(ag):

  with pytest.raises(Exception):
    ag.run(ngen=2,verbose=False,mode='standard',mutate_flip_pb=0.4,mutation_pb=1)

  ag.run(ngen=0,verbose=False,mode='standard',mutate_flip_pb=0.4,mutation_pb=1,crossover_pb=0)
  assert ag.ga.params['mutate_flip_pb'] == 0.4
  assert ag.ga.params['mutation_pb'] == 1
  assert ag.ga.params['crossover_pb'] == 0

def test_ga_individual_standard(ag):

  ag.run(ngen=0,verbose=False,mode='standard')

  for _ in range(100):
    ind = ag.ga.individual_standard()
    assert sum(ind) >= 5
    #print(ind.astype(int))
  
  ind = ag.ga.individual_standard()
  for _ in range(100):
    ind = ag.ga.mutate_standard(ind)[0]
    assert sum(ind) >= 5

  ind1 = ag.ga.individual_standard()
  ind2 = ag.ga.individual_standard()
  for _ in range(100):
    ag.ga.crossover_standard(ind1, ind2)
    assert sum(ind1) >= 5
    assert sum(ind2) >= 5

def test_ga_fixed(ag):

  ag.run(ngen=0,verbose=False,mode='fixed',nfeatures=5)
  for _ in range(100):
    ind = ag.ga.individual_fixed()
    assert sum(ind) == 5

  ind = ag.ga.individual_fixed()
  for _ in range(10):
    ind = ag.ga.mutate_fixed(ind)[0]
    #print(ind)
    assert sum(ind) == 5

  ind1 = ag.ga.individual_fixed()
  ind2 = ag.ga.individual_fixed()
  print(ind1.astype(int), ind2.astype(int))
  for _ in range(10):
    ag.ga.crossover_fixed(ind1, ind2)
    print(ind1.astype(int), ind2.astype(int))
    assert sum(ind1) == 5
    assert sum(ind2) == 5
    ag.ga.mutate_fixed(ind1)
    ag.ga.mutate_fixed(ind2)
    print(ind1.astype(int), ind2.astype(int))
    print()

def test_ga_helper_functions(ag):

  ag.run(ngen=0,verbose=False,mode='standard',mutate_flip_pb=0.4)

  #
  # Test fill_up
  #
  ind = np.full(10, False)
  ind[1] = True
  ag.ga.fill_up(ind)
  assert sum(ind) == 5

  #
  # Test n_bitflip
  #

  ind = ag.ga.individual_standard()
  for _ in range(10):
    old = ind.copy()
    print(old)
    ag.ga.n_bitflip(ind)
    print(ind)
    print(sum(ind != old))

  ind = ag.ga.individual_standard()
  ind.fill(False)
  ag.ga.gen.seed(0)
  ag.ga.n_bitflip(ind)
  assert np.array_equal(ind, [True, False, False, False, True, False, False, True, False, False])
  assert hasattr(ind, 'fitness')

  #
  # Test cross_and_or
  #
  ag.ga.gen.seed(0)
  ind1 = ag.ga.individual_standard()
  ind2 = ag.ga.individual_standard()

  ag.ga.cross_and_or(ind1,ind2)
  assert np.array_equal(ind1, [False, False, True, False, True, True, True, True, False, False])
  assert hasattr(ind1, 'fitness')

  #
  # Test swap_block
  #
  for _ in range(10):
    ind1 = np.full(10, True)
    ind2 = np.full(10, False)

    ag.ga.swap_block(ind1, ind2)

    assert sum(ind1) < 10
    assert sum(ind2) < 10
    assert sum(ind1) + sum(ind2) == 10
    assert np.array_equal(np.logical_or(ind1,ind2), np.full(10, True))
