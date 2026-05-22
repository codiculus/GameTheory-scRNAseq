# GameTheory-scRNAseq

This repository holds the code associated to the paper "Exploration of single-cell RNAseq expression variability considering gene cooperation", submitted to CDC 2026.

The project analizes the genes that contribute most to determinate cellular structures, introducing a new gene ranking in single-cell ribonucleic acid sequencing (scRNA-seq).

All the code is available in the src folder, which has been developed using Python which the library Cython for speeding up computations. Also, the scanpy library has been used to load GNA data.

Please refer to the requirements.txt file for a detailed list of dependencies.


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

### Get the sc-RNAseq data.

The original dataset is available from [Insert link here]()

### Run the algorithms

We provide the user with a script that generates the results and another that can be used to generate all the Figures presented in Section IV of the paper.

To run the classification and generate the Shapley values, please run:

```
 > ./generate_results.py
```

Once the results have been obtained, you can use the following jupyter Python notebooks:

```
 > ./violin_plots.ipynb
 > ./scatter_plots.ipynb
```
