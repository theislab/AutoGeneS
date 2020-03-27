import pytest

import sys
sys.path.insert(0,"..")

import numpy as np
from autogenes import AutoGeneS

@pytest.fixture
def ag_simple():
  data = np.zeros((3,6))
  data[0,0:2] = 1
  data[1,2:4] = 1
  data[2,4:6] = 1

  return AutoGeneS(data)

@pytest.fixture
def ag_random():
  np.random.seed(0)
  data = np.random.randint(0,5,(3,10))
  return AutoGeneS(data)

def test_run(ag_simple):

  ag_simple.run(ngen=10,seed=0,mode='fixed',nfeatures=3,verbose=False,weights=(1,-1),objectives=('distance','correlation'))
  
  s = ag_simple.select(weights=(1,0))
  assert np.array_equal([False, True, True, False, False, True], s)

def test_from_pareto(ag_random):

  ag_random.run(ngen=10,verbose=False,weights=(1,-1),objectives=('distance','correlation'))

  fit = ag_random.fitness_matrix
  print(fit)
  
  # The selection by weight uses the normalized fitness matrix
  fit_n = fit.copy()
  for i in range(fit_n.shape[1]):
    fit_n[:,i] *= 1/np.max(fit_n[:,i])
  print(fit_n)

  with pytest.raises(Exception, match="You need to provide exactly one criterion."):
    ag_random._AutoGeneS__from_pareto()

  with pytest.raises(Exception, match="You need to provide exactly one criterion."):
    ag_random._AutoGeneS__from_pareto(weights=(1,1),close_to=12)

  assert ag_random._AutoGeneS__from_pareto(weights=(1,0)) == (0,'weights')
  assert ag_random._AutoGeneS__from_pareto(weights=(-1,0)) == (4,'weights')
  assert ag_random._AutoGeneS__from_pareto(weights=(1,-1)) == (2,'weights')

  assert ag_random._AutoGeneS__from_pareto(index=2) == (2,'index')
  with pytest.raises(ValueError, match="Invalid index"):
    ag_random._AutoGeneS__from_pareto(index=10)
  assert ag_random._AutoGeneS__from_pareto(index=(0,0)) == (4,'index')
  assert ag_random._AutoGeneS__from_pareto(index=(0,1)) == (3,'index')
  assert ag_random._AutoGeneS__from_pareto(index=(1,4)) == (0,'index')
  assert ag_random._AutoGeneS__from_pareto(index=(1,-1)) == (0,'index')

  assert ag_random._AutoGeneS__from_pareto(close_to=(0,17)) == (0,'close_to')
  assert ag_random._AutoGeneS__from_pareto(close_to=(1,0.32)) == (3,'close_to')
  assert ag_random._AutoGeneS__from_pareto(close_to=(1,0)) == (4,'close_to')

def test_select(ag_random):

  ag = ag_random

  ag.run(ngen=10,verbose=False,weights=(1,-1),objectives=('distance','correlation'))
  assert ag.select(close_to=(0,17)) == ag.pareto[0]
  assert ag.selection == ag.pareto[0]
  assert ag.selection_index == 0

  assert ag.select(index=(0,1)) == ag.pareto[3]
  assert ag.selection == ag.pareto[3]
  assert ag.selection_index == 3
