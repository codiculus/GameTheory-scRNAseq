# GameTheory-scRNAseq

This repository holds the code associated to the paper "Exploration of single-cell RNAseq expression variability considering gene cooperation", submitted to ECC 2026.

The project analizes the genes that contribute most to determinate cellular structures.

All the code is available in the src folder, which has been developed using Python which the library Cython for speeding up computations. Also, the scanpy library has been used to load GNA data.


## Installation

## Usage


### Dependencies

Our code runs in Python 3.12 or greater, with the libraries specified in the requirements.txt file.

### Compilation

Go to the src subfolder and compile the Cython code to obtain the executable files:

```
 > cd src
 > python setup.py build_ext --inplace
```

### Get the GNA seq data.


### Run the algorithms

Run the Jupyter notebook to obtain the Figures  .... of the paper.


