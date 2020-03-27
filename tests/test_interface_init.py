import pytest

import sys
sys.path.insert(0,"..")

import autogenes as ag
import numpy as np
import pandas as pd
import anndata

def test_compute_means():

  ann = anndata.AnnData(np.reshape(np.arange(70),(7,10)))
  ann.obs['ct'] = ['H1','H2','H1','G','H1','H2','H2']

  result = ag.main._Interface__compute_means(ann,'ct')
  assert result.X.shape == (3,10)
  assert result.var_vector('H1')[0] == 20

def test_interface_init():

  with pytest.raises(Exception):
    ag.init(None)

  # numpy arrays
  data = np.identity(3)
  ag.init(data)
  assert ag.main.data is data
  assert ag.main.main.data is data
  assert sum(ag.main.pre_selection) == 3

  with pytest.raises(Exception):
    ag.init(np.zeros((5,2)))

  # DataFrame
  cols = ["gene_1","gene_2","gene_3"]
  df = pd.DataFrame(data,columns=cols)
  ag.init(df)

  assert np.all(ag.main.data == df.values)
  assert np.all(ag.main.main.data == df.values)
  assert np.all(ag.main.data_genes == df.columns.values)
  assert sum(ag.main.pre_selection) == 3

  # Check reset
  ag.init(data)
  assert ag.main.data_genes is None

  # AnnData
  test_data = np.zeros((7,5))
  test_data[[0,2,4],:] = 1
  test_data[[1,5,6],:] = -3
  test_data[3,:] = 4
  adata = anndata.AnnData(test_data)
  genes = [f"gene_{i}" for i in range(5)]
  adata.var_names = genes
  adata.var["highly_variable"] = np.array([True,True,False,True,False])
  adata.var["selection"] = np.array([False,False,True,True,True])
  adata.obs["col"] = np.full((7,), True)

  # No celltype column
  with pytest.raises(ValueError):
    ag.init(adata)

  adata.obs["celltype"] = ['H1','H2','H1','G','H1','H2','H2']
  adata.obs["only_h1"] = ['H1','H1','H1','H1','H1','H1','H1']
  adata.obs["celltype2"] = ['H1','H2','H1','H1','H1','H2','H2']
  
  # Simple
  ag.init(adata)

  test_data_mean = np.repeat(np.array([4,1,-3]).reshape(3,1), 5, axis=1)

  assert np.array_equal(ag.main.data_genes, genes)
  assert sum(ag.main.pre_selection) == 5
  assert np.array_equal(ag.main._adata.var_names, genes)
  assert np.array_equal(ag.main._adata.obs_names, ['G', 'H1','H2'])
  assert np.array_equal(ag.main._adata.X, test_data_mean)
  assert np.array_equal(ag.main.data, test_data_mean)

  # celltype_key
  with pytest.raises(ValueError):
    ag.init(adata,celltype_key="only_h1")

  ag.init(adata,celltype_key="celltype2")
  assert np.array_equal(ag.main._adata.X, np.repeat(np.array([1.75,-3]).reshape(2,1), 5, axis=1))

  # genes_key
  ag.init(adata,genes_key="selection")
  assert np.array_equal(ag.main._adata.X, test_data_mean)
  assert np.array_equal(ag.main.data, test_data_mean[:,[2,3,4]])
  assert np.all(ag.main.pre_selection == adata.var["selection"])
  # Not the selected genes, but ALL original genes!
  assert np.all(ag.main.data_genes == adata.var_names.values)

  # use_highly_variable
  ag.init(adata,use_highly_variable=True)
  assert np.array_equal(ag.main._adata.X, test_data_mean)
  assert np.array_equal(ag.main.data, test_data_mean[:,[0,1,3]])
  assert np.all(ag.main.pre_selection == adata.var["highly_variable"])
  assert np.all(ag.main.data_genes == adata.var_names.values)

  # celltype_key (2)
  test_data = np.reshape(np.arange(70),(7,10))
  ann = anndata.AnnData(test_data)
  ann.obs["ct"] = ['H1','H2','H1','G','H1','H2','H2']

  ag.init(ann,celltype_key="ct")
  assert np.all(ag.main._adata.shape == (3,10))
  assert np.all(ag.main._adata.var_vector('H1') == [20+i for i in range(10)])
  assert np.all(ag.main._adata.var_vector('G') == [30+i for i in range(10)])
  assert np.all(ag.main.data_genes == ann.var_names.values)
  assert sum(ag.main.pre_selection) == 10

  # celltype_key + genes_key
  sel = np.array([True,True,True,False,False,True,False,True,True,True])
  ann.var['selection'] = sel
  ag.init(ann,celltype_key="ct",genes_key="selection")

  # NOT (3,7)! The genes are not applied, but stored in pre_selection
  assert ag.main._adata.X.shape == (3,10)
  assert ag.main.data.shape == (3,7)

  sel_ids, = np.where(sel)
  assert np.array_equal(ag.main.data[1], [20+i for i in sel_ids])
  assert np.array_equal(ag.main.data[0], [30+i for i in sel_ids])
  assert np.array_equal(ag.main.data_genes,ann.var_names)
  assert sum(ag.main.pre_selection) == 7
