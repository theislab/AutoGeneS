import pytest

import numpy as np
import pandas as pd

from scipy.special import binom

import os
import sys
sys.path.insert(0, "..")

from autogenes import objectives as ga_objectives

def test_distance():

  arr = np.ones((3,3))
  assert ga_objectives.distance(arr) == 0

  arr = np.identity(3)
  assert np.isclose(ga_objectives.distance(arr), 3 * np.sqrt(2))

def test_correlation():

  arr = np.ones((3,3))
  # Should'nt throw a warning
  assert ga_objectives.correlation(arr) == 0

  arr = np.zeros((3,3))
  assert ga_objectives.correlation(arr) == 0

  arr = np.identity(5)
  # Let i!=j
  # cov(e_i, e_j) = 1/4 * ( 2 * (4/5) * (-1/5) + 3 * (-1/5) * (-1/5)) = -1/20
  # var(e_i) = 1/4 * ( (4/5)^2 + 4 * (-1/5)^2) = 1/5
  # (i,j)th entry in the corrcoef matrix is (-1/20)/sqrt(var(e_i) * var(e_j)) = -1/4
  # Result is 1/4 * (5 over 2)

  assert np.isclose(ga_objectives.correlation(arr), 1/4*binom(5,2))
