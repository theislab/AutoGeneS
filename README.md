# AutoGeneS

Automatic Gene Selection (for bulk deconvolution)

## Dependencies

* python >=3.6
* dill
* deap
* cachetools

## Installation
git clone https://github.com/theislab/AutoGeneS<br/>
pip install --user dist/autogenes-0.9.1-py3-none-any.whl<br/>

## Testing
import numpy as np<br/>
from autogenes import AutoGenes<br/>
ag = AutoGenes(np.identity(2))<br/>

## Examples

Simple examples for how to call AutoGeneS can be find in folder "jupyter"

Bulk deconvolution of synthetic bulks using single-cell reference profiles can be find in folder "deconv_example"
