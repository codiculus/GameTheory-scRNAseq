import numpy as np
import itertools
import math

def calcular_CV_grupo_genes(m, posiciones_genes):
    
    # Inicializar la suma con ceros
    suma = np.zeros(m.shape[1])

    for gen in posiciones_genes:
        
        valores_gen = m[gen, :]  # Obtener los valores de expresión del gen
        suma += valores_gen

    varianza = np.var(suma)
    media = np.mean(suma)

    coeficiente_variacion = varianza/media

    return coeficiente_variacion

def valores_shapley(jugadores, valores_coaliciones):
    """
    Calcula los valores de Shapley para los jugadores.
    
    Args:
        players (list): Lista de jugadores (ej: ["A","B","C"]).
        coalition_values (dict): Diccionario con los valores de las coaliciones.
                                 La clave es un frozenset de jugadores, ej:
                                 {frozenset(["A"]): 10, frozenset(["B","C"]): 25, ...}

    Returns:
        dict: Valor de Shapley de cada jugador.
    """
    n = len(jugadores)
    shapley = {p: 0 for p in jugadores}

    for i in jugadores:
        # recorrer todos los subconjuntos posibles que no incluyan a i
        for r in range(n):
            for S in itertools.combinations([p for p in jugadores if p != i], r):
                S = set(S)
                S_frozen = frozenset(S)
                S_with_i = frozenset(S | {i})

                v_S = valores_coaliciones.get(S_frozen, 0)
                v_S_i = valores_coaliciones.get(S_with_i, 0)

                # peso combinatorio
                weight = math.factorial(len(S)) * math.factorial(n - len(S) - 1) / math.factorial(n)

                shapley[i] += weight * (v_S_i - v_S)

    return shapley