from math import exp
import numpy as np

def ganancias_cython(combinaciones, float[:,:] ganancias):

    cdef int N = len(combinaciones)
    cdef double[:] Y = np.zeros(N)
    cdef int i, j, d
    cdef double r = 0

    for i in range(N):
        print (combinaciones[i])
        Y[i] = calcular_valores_coalicion_cython(ganancias, combinaciones[i])

    return Y

cdef calcular_valores_coalicion_cython(float[:,:] m, int [:]pos_genes):
    """
    Calcula el valor de una coalición de genes usando diferentes métodos.
    
    Args:
        m (numpy.array): Matriz de expresión (genes x células)
        pos_genes (list): Índices de los genes a sumar
        normalizar (bool): aplica normalización después de sumar.
                          
        metodo (str): Método de cálculo.
                     'varianza' - Varianza simple
    
    Returns:
        float: Valor calculado según el método seleccionado
    """
    # Inicializar la suma con ceros
    cdef double [:]suma = np.zeros(m.shape[1])
    cdef int i
    cdef int [:]total_counts_per_cell
    cdef double resultado

    # Sumar la expresión de los genes especificados
    for i in range(len(pos_genes)):
        valores_gen = m[pos_genes[i], :]  # Obtener los valores de expresión del gen
        for j in range(m.shape[1]):
            suma[j] += valores_gen[j]

    
    # NORMALIZACIÓN: Aplicar normalización después de la suma
    # 1. Normalizar por número total de counts por célula
    total_counts_per_cell = np.sum(m, axis=0)  # Suma total por célula
    for i in range(m.shape[1]):
        suma[i] = 1 + (suma[i] / total_counts_per_cell[i]) * 10000  # target_sum = 1e4
    
    
    # 2. Aplicar transformación logarítmica
    suma = np.log2(suma)  # log(x + 1)

    # Calcular según el método seleccionado
    resultado = np.var(suma)
    
    return resultado