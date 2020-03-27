from .interface import Interface
from .core import AutoGeneS

main = Interface()

init = main.init
optimize = main.optimize
select = main.select
plot = main.plot
deconvolve = main.deconvolve
pipeline = main.pipeline
pareto = main.pareto
resume = main.resume
save = main.save
load = main.load
adata = main.adata
fitness_matrix = main.fitness_matrix
selection = main.selection
