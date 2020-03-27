import pytest

import sys
sys.path.insert(0,"..")

import autogenes as ag
import numpy as np
import pandas as pd
import anndata

def test_pareto_fitness():

  ag.init(np.zeros((3,10)))
  ag.optimize(mode='fixed',nfeatures=5,verbose=False,weights=(1,-1),objectives=('distance','correlation'))
  assert np.array_equal(ag.pareto(), ag.main.main.pareto)

  np.random.seed(20)
  adata = anndata.AnnData(np.random.randint(-5,5,(2,10)))
  adata.obs['celltype'] = ['1','2']
  adata.var['highly_variable'] = np.array([True,True,True,True] + [False]*6)
  ag.init(adata,use_highly_variable=True)
  ag.optimize(ngen=3,mode='fixed',nfeatures=2,verbose=False,weights=(1,-1),objectives=('distance','correlation'))
  for p in ag.pareto():
    assert np.array_equal(p[4:], [False]*6)

  print(ag.pareto())
  print(ag.fitness_matrix())

  assert ag.fitness_matrix().shape == (len(ag.pareto()), 2)

def test_process_selection():

  ag.main.pre_selection = np.array([True, True, False, True, False, True, True])

  process_selection = ag.main._Interface__process_selection
  
  s1 = np.array([False,True,True,False,False])
  s2 = np.array([True,True,False,False,False])
  s3 = np.array([True,True,False,False,False,False])
  
  assert np.all(process_selection(s1) == [False,True,False,True,False,False,False])
  assert np.all(process_selection(s2) == [True,True,False,False,False,False,False])
  with pytest.raises(Exception):
    assert process_selection(s3)

def test_select():

  data = np.zeros((3,6))
  data[0,0:2] = 1
  data[1,2:4] = 1
  data[2,4:6] = 1

  adata = anndata.AnnData(data)
  adata.obs['celltype'] = ['c1','c2','c3']
  adata.var['highly_variable'] = [True,True,True,True,False,False]

  ag.init(adata,use_highly_variable=True)
  ag.optimize(ngen=10,seed=0,mode='fixed',nfeatures=3,verbose=False,weights=(1,-1),objectives=('distance','correlation'))
  
  s = ag.main.main.select(weights=(1,0))
  s_target = [True, True, False, True]
  s_target_processed = [True,True,False,True,False,False]
  assert np.array_equal(s_target, s)
  assert np.array_equal(s_target_processed, ag.main._Interface__process_selection(s))

  res = ag.select(weights=(1,0))
  assert np.array_equal(ag.adata().var['autogenes'], s_target_processed)
  assert np.array_equal(res, s_target_processed)
  assert np.array_equal(ag.selection(), s_target_processed)

  res_adata = ag.select(weights=(1,0), copy=True)
  assert not (res_adata is ag.main.adata)
  assert np.array_equal(res_adata.var['autogenes'], s_target_processed)

def test_select2():

  np.random.seed(0)
  adata = anndata.AnnData(np.random.randint(-5,5,(3,5)))
  adata.obs['celltype'] = ['c1','c2','c3']
  
  ag.init(adata)
  ag.optimize(ngen=10,offspring_size=100,verbose=False,mode='fixed',nfeatures=3,weights=(1,-1),objectives=('distance','correlation'))

  ag.select(weights=(1,-1),key_added='autogenes',copy=False)
  ag.select(weights=(1,0),key_added='autogenes2',copy=False)
  assert not np.array_equal(ag.adata().var['autogenes'], ag.adata().var['autogenes2'])

  r1 = ag.select(weights=(1,-1),key_added='autogenes',copy=True)
  r2 = ag.select(weights=(1,0),key_added='autogenes2',copy=True)

  assert not ('autogenes2' in r1.var_names)
  assert not ('autogenes' in r2.var_names)
