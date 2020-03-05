.. AutoGenes documentation master file, created by
   sphinx-quickstart on Tue Mar  3 22:08:44 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AutoGenes
==========

AutoGeneS automatically extracts informative genes and reveals the cellular heterogeneity of bulk RNA samples. AutoGeneS requires no prior knowledge about marker genes and selects genes by simultaneously optimizing multiple criteria: minimizing the correlation and maximizing the distance between cell types. It can be applied to reference profiles from various sources like single-cell experiments or sorted cell populations.

For a multi-objective optimization problem, there usually exists no single solution that simultaneously optimizes all objectives. In this case, the objective functions are said to be conflicting, and there exists a (possibly infinite) number of Pareto-optimal solutions. Pareto-(semi)optimal solutions are a set of all solutions that are not dominated by any other explored solution. Pareto-optimal solutions offer a set of equally good solutions from which to select, depending on the dataset

.. toctree:: 
   :caption: Contents:
   :maxdepth: 2
   installation
   tutorial
