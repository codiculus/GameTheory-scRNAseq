import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import itertools
import gametheory as gt
import csv
import time

inicio = time.time()

ruta_datos = "../filtered_gene_bc_matrices/hg19/"

adata = sc.read_10x_mtx(ruta_datos, var_names='gene_symbols')

print(adata)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(adata)


# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)



adata = adata[
    (adata.obs.n_genes_by_counts < 2500)
    & (adata.obs.n_genes_by_counts > 200)
    & (adata.obs.pct_counts_mt < 5),
    :,
].copy()
adata.layers["counts"] = adata.X.copy()


sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

print(adata)

sc.pp.highly_variable_genes(adata,layer="counts",n_top_genes=20,min_mean=0.0125,max_mean=3,min_disp=0.5,flavor="seurat_v3")
sc.pl.highly_variable_genes(adata)

print(adata.var['highly_variable'])

genes_mas_variables = adata.var[adata.var['highly_variable']]

indices_genes_mas_variables = np.where(adata.var['highly_variable'])[0]

#Guardar genes mas variables
genes_mas_variables.to_csv("./results/genes_mas_variables.csv")
np.savetxt("./results/indices_genes_mas_variables.csv", indices_genes_mas_variables, delimiter=",")

genes_mas_variables_array = np.array(genes_mas_variables)

nombres_genes_mas_variables = adata.var[adata.var['highly_variable']].index.tolist()

matriz_norm, valores_coaliciones = gt.calcula_valores_coaliciones_cython(adata, nombres_genes_mas_variables, indices_genes_mas_variables)

#print(valores_coaliciones)

y = gt.valores_shapley(indices_genes_mas_variables, {frozenset(k): v for k, v in valores_coaliciones.items()})

print("Shapley Values:")

for jugador, valor_shap in y.items(): 
    print(f"Player: {jugador}, Shapley Value: {valor_shap}")

fin = time.time()
tiempo_total = fin - inicio

horas = int(tiempo_total // 3600)
minutos = int((tiempo_total % 3600) // 60)
segundos = tiempo_total % 60

print(f"Execution time: {horas}h {minutos}m {segundos:.2f}s")
