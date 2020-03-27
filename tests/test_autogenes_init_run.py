import pytest

import sys
sys.path.insert(0,"..")

import numpy as np
from autogenes import AutoGeneS
from autogenes import objectives as ga_objectives

def test_autogenes_init():
  
  with pytest.raises(Exception):
    AutoGeneS()

  with pytest.raises(Exception):
    AutoGeneS([[11,1],[1]])

  with pytest.raises(ValueError, match="data is expected to have two dimensions"):
    AutoGeneS(np.array([1]))

  with pytest.raises(ValueError, match="At least two rows \\(cell types\\) expected"):
    AutoGeneS(np.int_([[1,2,3]]))

  with pytest.raises(ValueError, match="Number of columns \\(genes\\) must be >= number of rows \\(cell types\\)"):
    AutoGeneS(np.int_(np.zeros((5,2))))

  with pytest.raises(ValueError, match="Number of columns \\(genes\\) must be >= number of rows \\(cell types\\)"):
    AutoGeneS(np.int_(np.zeros((2,1))))

  with pytest.raises(ValueError, match="Some entries in data are not scalars"):
    AutoGeneS(np.array([[1,np.inf],[1,np.nan]]))

  arr1 = np.float64([[-1, 12.5], [1, 1E10]])
  ag1 = AutoGeneS(arr1)
  assert ag1.selection == None
  assert np.array_equal(ag1.data,arr1)

def test_autogenes_run_general_args():
  
  ag = AutoGeneS(np.random.randn(2,5))

  # mode = standard

  ag.run(ngen=0,mode='standard',verbose=False)

  # mode = fixed

  with pytest.raises(ValueError, match="You need to supply nfeatures"):
    ag.run(mode='fixed')

  with pytest.raises(ValueError, match="nfeatures doesn't apply to standard mode \\(did you mean mode='fixed'\\?\\)"):
    ag.run(nfeatures=10)

  with pytest.raises(ValueError, match="nfeatures must be <= the number of columns \\(genes\\)"):
    ag.run(mode='fixed', nfeatures=10)

  with pytest.raises(ValueError, match="nfeatures must be >= the number of rows \\(cell types\\)"):
    ag.run(mode='fixed', nfeatures=1)

  ag.run(ngen=0, mode='fixed', nfeatures=2, verbose=False)

  with pytest.raises(ValueError, match="No such objective"):
    ag.run(mode='fixed', nfeatures=3, weights=(1,1), objectives=('a','b'))

def test_autogenes_run_objectives_and_weights():

  ag = AutoGeneS(np.random.randn(2,5))

  ag.run()
  assert ag.objectives_func == [ga_objectives.correlation, ga_objectives.distance]

  with pytest.raises(Exception, match="Need objectives for weights"):
    ag.run(weights=(1,1))

  with pytest.raises(Exception, match="Need weights for objectives"):
    ag.run(objectives=('distance','test'))

  with pytest.raises(ValueError, match="No such objective: test"):
    ag.run(weights=(1,1), objectives=('distance','test'))

  with pytest.raises(ValueError, match="Invalid objective"):
    ag.run(weights=(1,1), objectives=(1,1))

  with pytest.raises(ValueError, match="Number of weights does not match number of objectives"):
    ag.run(weights=(3,4,5),objectives=('distance','correlation'))

  try:
    ag.run(weights=(1,),objectives=('distance',),verbose=False)
  except Exception:
    pytest.fail("Should support single objective")

  with pytest.raises(Exception):
    ag.run(weights=1, objectives=('distance',))

  with pytest.raises(Exception):
    ag.run(weights=('a',),objectives=('distance',))

  with pytest.warns(UserWarning):
    ag.run(weights=(3,0,5),objectives=('distance','num_genes','correlation'))
    assert ag.weights == (3,5)
    assert ag.objectives_func == [ga_objectives.distance,ga_objectives.correlation]
    assert ag.objectives_names == ['distance','correlation']

def test_autogenes_run_custom_objectives():

  ag = AutoGeneS(np.random.randn(2,5))

  # Lambda as objective
  obj1 = lambda x: len(x)
  ag.run(ngen=0,weights=(1,1), objectives=(obj1, 'distance'), verbose=False)

  assert ag.objectives_func == [obj1, ga_objectives.distance]
  assert ag.objectives_names == ['<lambda>', 'distance']

  # Invalid objective
  def invalid_objective(data): 
    return data

  with pytest.raises(Exception):
    ag.run(ngen=0,weights=(-1,1), objectives=(invalid_objective, 'distance'), verbose=False)

  # Normal function
  def count(data): 
    return data.shape[0]

  ag.run(ngen=0,weights=(-1,1), objectives=(count, 'distance'), verbose=False)
  assert ag.objectives_func == [count, ga_objectives.distance]
  assert ag.objectives_names == ['count', 'distance']
