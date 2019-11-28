# AutoGeneS

Automatic Gene Selection (for bulk deconvolution)

## Dependencies

* python >=3.6
* anndata
* deap
* cachetools

## Installation
git clone https://github.com/theislab/AutoGeneS<br/>
pip install --user package/v0.9.1/dist/autogenes-0.9.1-py3-none-any.whl<br/>

## Testing
import numpy as np<br/>
from autogenes import AutoGenes<br/>
ag = AutoGenes(np.identity(2))<br/>

## Examples

Find examples here: "AutoGeneS/v0.9.1/jupyter"