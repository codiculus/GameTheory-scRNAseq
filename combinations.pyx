import numpy as np
from array import array
from cython.parallel import prange
from cython.cimports.libc.math import log2
import cython
from libc.stdlib cimport malloc, free

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.

def ganancias_cython(combinaciones, float[:,:] m, int max_cols):

    cdef int N = len(combinaciones)
    cdef int M = m.shape[1]
    cdef double[:] Y = np.zeros(N)
    cdef int i, j, k, l
    cdef double r = 0
    cdef double mean
    cdef double [:]suma = np.zeros(M)
    cdef double var
    cdef float [:]total_counts_per_cell = np.sum(m, axis=0)  # Suma total por célula
    cdef int c_array[20000]
    cdef int c_matrix[1024][10]
    #cdef int *c_array = <int *> malloc(N * max_cols * sizeof(int))    

    for i in range(N):
        limite = len(combinaciones[i])
        for j in range(max_cols):
            if j < limite:
                c_array[i*N + j] = combinaciones[i][j]
                c_matrix[i][j] = combinaciones[i][j]
                
            else:
                c_matrix[i][j] = -1
                c_array[i*N + j] = -1

            print (c_matrix[i][j], c_array[i*N+j])

        

    

    print("Combinaciones calculadas")
    for i in prange(N, nogil=True):
        for l in range(M):
            suma[l] = 0.0
        
        # Sumar la expresión de los genes especificados
        for j in range(max_cols):
            k = c_array[i*N + j]
            k = c_matrix[i][j]
            if k < 0 or k >= 13714:
                break
            for l in range(M):
                suma[l] += m[k, l]


        mean = 0.0
        for l in range(M):
            suma[l] = log2(1 + (suma[l] / total_counts_per_cell[l]) * 10000)
            mean = mean + suma[l]
        mean = mean / M
        var = 0.0
        for l in range(M):
            var = var + (suma[l] - mean) * (suma[l] - mean)
        Y[i] = var / M
    
    return Y


