import scipy.spatial.distance
import numpy as np

def distance(data):

  total = 0
  num_obs = data.shape[0]

  for i in range(num_obs):
    for j in range(i+1,num_obs):
      total += scipy.spatial.distance.euclidean(data[i],data[j])

  return total

def correlation(data):

  total = 0 

  with np.errstate(divide='ignore',invalid='ignore'):
    corr = np.corrcoef(data)

  # xij = NaN only occurs if either the variable xi or xj had variance 0
  # Then, cov(xi, xj) = 0
  # Therefore, we can set NaN values to 0
  corr[np.isnan(corr)] = 0

  num_vars = corr.shape[0]
  for i in range(num_vars):
    for j in range(i+1,num_vars):
      total += abs(corr[i,j])

  return total

def condition(data):
  return np.linalg.cond(data, p=2)
