import pytest

import sys
sys.path.insert(0,"..")

import autogenes as ag
import numpy as np
import pandas as pd
import anndata

from sklearn.svm import NuSVR
from sklearn import linear_model
from scipy.optimize import nnls

def test_unpack_bulk():

  unpack_bulk = ag.main._Interface__unpack_bulk

  arr = np.ones((3,))
  bulk_data, bulk_genes = unpack_bulk(arr)
  assert np.array_equal(bulk_data,arr.reshape(1,3))
  assert type(bulk_data) == np.ndarray

  arr = np.ones((2,3))
  bulk_data, bulk_genes = unpack_bulk(arr)
  assert bulk_data is arr
  assert bulk_genes is None
  assert type(bulk_data) == np.ndarray

  adata = anndata.AnnData(arr)
  gene_names =["gene_1","gene_2","gene_3"]
  adata.var_names = gene_names
  bulk_data, bulk_genes = unpack_bulk(adata)
  assert np.array_equal(bulk_data,adata.X)
  assert np.array_equal(bulk_genes, gene_names)
  assert type(bulk_data) == np.ndarray
  assert type(bulk_genes) == np.ndarray

  series = pd.Series([1,2,3],index=gene_names)
  bulk_data, bulk_genes = unpack_bulk(series)
  assert np.array_equal(bulk_data,series.values.reshape(1,3))
  assert np.array_equal(bulk_genes, gene_names)
  assert type(bulk_data) == np.ndarray
  assert type(bulk_genes) == np.ndarray

  df = adata.to_df()
  bulk_data, bulk_genes = unpack_bulk(df)
  assert np.array_equal(bulk_data,df.values)
  assert np.array_equal(bulk_genes, gene_names)
  assert type(bulk_data) == np.ndarray
  assert type(bulk_genes) == np.ndarray

def test_model_input():

  model_input = ag.main._Interface__model_input

  # Note: all bulk_data must be 2-dim since this is the output of __unpack_bulk

  #
  # No gene labels
  #
  
  # Shape mismatch
  ag.main.data = np.zeros((2,3))
  ag.main.data_genes = None
  bulk_data = np.zeros((2,2))
  bulk_genes = None
  with pytest.raises(ValueError):
    model_input(bulk_data,bulk_genes,np.full((3,),True))

  # Simple case (with 1-dim bulk_data)
  ag.main.data = np.reshape(np.arange(6),(2,3))
  ag.main.data_genes = None
  ag.main.selection = np.full((3,),True)
  bulk_data = np.array([1,2,3]).reshape(1,3)
  bulk_genes = None

  X,y = model_input(bulk_data,bulk_genes,ag.main.selection)
  assert np.array_equal(X,ag.main.data.T)
  assert np.all(y == bulk_data.reshape(3,1))

  # With selection
  ag.main.selection = np.array([True,True,False])
  X,y = model_input(bulk_data,bulk_genes,ag.main.selection)
  assert np.array_equal(X,ag.main.data.T[0:2])
  assert np.all(y == bulk_data[:,0:2].reshape(2,1))

  #
  # With gene labels
  #

  gene_names = np.array(["gene_1","gene_2","gene_3"])

  # Simple case
  ag.main.data = np.reshape(np.arange(6),(2,3))
  ag.main.data_genes = gene_names
  ag.main.selection = np.array([True,True,False])
  bulk_data = np.reshape(np.arange(6),(2,3))
  bulk_genes = gene_names

  X,y = model_input(bulk_data,bulk_genes,ag.main.selection)
  assert np.array_equal(X,ag.main.data.T[0:2])
  assert np.all(y == bulk_data.T[0:2])

  # Permutation
  bulk_genes = np.array(["gene_2","gene_1","gene_3"])
  X,y = model_input(bulk_data,bulk_genes,ag.main.selection)
  assert np.array_equal(X,ag.main.data.T[0:2])
  assert np.all(y == bulk_data.T[[1,0]])

  # Complex case
  ag.main.data = np.reshape(np.arange(8),(2,4))
  ag.main.data_genes = np.array(["gene_1","gene_2","gene_3","gene_4"])
  ag.main.selection = np.array([True,True,True,False])
  bulk_data = np.reshape(np.arange(4),(2,2))
  bulk_genes = np.array(["gene_3","gene_2"])

  with pytest.warns(UserWarning):
    X,y = model_input(bulk_data,bulk_genes,ag.main.selection)
  assert np.array_equal(X,ag.main.data.T[1:3])
  assert np.all(y.T == bulk_data[:,[1,0]])

def test_deconvolve():

  data = pd.read_csv('../datasets/GSE75748_bulk_data.csv',index_col='index')
  data = data.T.iloc[:,:10]
  data.head()

  # Basic tests of models
  ag.init(data)
  ag.optimize(seed=0)
  ag.select()

  r1 = ag.deconvolve(data.values[3],model='linear')
  idm = np.identity(6)
  assert np.allclose(r1[0], idm[3])

  X = data.T[ag.selection()]
  y = data.values[3,ag.selection()]

  r2 = ag.deconvolve(data.values[3],model='nusvr')
  model = NuSVR(nu=0.5,C=0.5,kernel='linear')
  model.fit(X, y)
  assert np.array_equal(model.coef_,r2)

  r3 = ag.deconvolve(data.values[3],model='nnls')
  x,err = nnls(X,y)
  assert np.array_equal(r3, [x])

  # Multiple targets
  A = np.array([[1,1,2],[2,0,1],[0,3,0]])
  bulk = A.dot(data.values[:3])
  r1 = ag.deconvolve(bulk,model='linear')
  target = np.hstack((A,np.zeros((3,3))))
  assert np.allclose(r1,target)

  #
  # Different bulk data types
  #

  adata = anndata.AnnData(data)
  adata.obs['celltype'] = [str(i) for i in range(data.shape[0])]
  ag.init(adata)
  ag.optimize(ngen=3)
  s = ag.select(index=3)
  print("All genes:", data.columns)
  print("Selected genes:", data.columns.values[s])

  # AnnData + AnnData
  bulk = adata[1:3]
  print(ag.deconvolve(bulk,model='linear'), np.identity(data.shape[0])[1:3])
  assert np.allclose(ag.deconvolve(bulk,model='linear'), np.identity(data.shape[0])[1:3], atol=1E-2)

  # AnnData + Series
  bulk_series = pd.Series(data=[1,1,1,2],index=['APOA2', 'OR9A4', 'TMIE','PHYHIPL' ])
  with pytest.warns(UserWarning):
    ag.deconvolve(bulk_series,model='linear')

  # Note: Force a minimal number of intersecting genes? (e.g. >= number of cell types?)
  bulk_series = pd.Series(data=[1],index=['TMIE' ])
  with pytest.warns(UserWarning):
    ag.deconvolve(bulk_series,model='linear')

  bulk_series = pd.Series(data=[1],index=['TMIE2' ])
  with pytest.raises(ValueError):
    ag.deconvolve(bulk_series,model='linear')

  bulk_series = data[['APOA2','OR9A4', 'TMIE','PHYHIPL']].iloc[1,:]
  with pytest.warns(UserWarning):
    assert np.allclose(ag.deconvolve(bulk_series,model='linear'), np.identity(data.shape[0])[1], atol=1E-2)

  # AnnData + DataFrame
  bulk_df = data.iloc[:2,:].copy()
  bulk_df.iloc[1] = bulk_df.iloc[0] * 2 + bulk_df.iloc[1]

  target = np.identity(data.shape[0])[:2]
  target[1] = target[0] * 2 + target[1]

  assert np.allclose(ag.deconvolve(bulk_df,model='linear'), target, atol=1E-2)
