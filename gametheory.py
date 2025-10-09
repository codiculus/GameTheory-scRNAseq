import numpy as np
import itertools
import math

def calcular_valores_coaliciones(m, posiciones_genes, normalizar=True, metodo='cv'):
    """
    Calcula el valor de una coalición de genes usando diferentes métodos.
    
    Args:
        m (numpy.array): Matriz de expresión (genes x células)
        posiciones_genes (list): Índices de los genes a sumar
        normalizar (bool): Si True, aplica normalización después de sumar.
                          Si False, usa los datos tal como están.
        metodo (str): Método de cálculo. Opciones:
                     'cv' - Coeficiente de variación (varianza/media)
                     'varianza' - Varianza simple
    
    Returns:
        float: Valor calculado según el método seleccionado
    """
    # Inicializar la suma con ceros
    suma = np.zeros(m.shape[1])

    # Sumar la expresión de los genes especificados
    for gen in posiciones_genes:
        valores_gen = m[gen, :]  # Obtener los valores de expresión del gen
        suma += valores_gen

    if normalizar:
        # NORMALIZACIÓN: Aplicar normalización después de la suma
        # 1. Normalizar por número total de counts por célula
        total_counts_per_cell = np.sum(m, axis=0)  # Suma total por célula
        suma_procesada = (suma / total_counts_per_cell) * 10000  # target_sum = 1e4
        
        # 2. Aplicar transformación logarítmica
        suma_procesada = np.log1p(suma_procesada)  # log(x + 1)
    else:
        # SIN NORMALIZACIÓN: Usar los datos tal como están
        suma_procesada = suma
    
    # Calcular según el método seleccionado
    if metodo == 'cv':
        # Calcular coeficiente de variación
        varianza = np.var(suma_procesada)
        media = np.mean(suma_procesada)
        resultado = varianza/media if media > 0 else 0
    elif metodo == 'varianza':
        # Calcular varianza simple
        resultado = np.var(suma_procesada)
    else:
        raise ValueError(f"Método '{metodo}' no reconocido. Use 'cv' o 'varianza'.")

    return resultado

def valores_shapley(jugadores, valores_coaliciones):
    """
    Calcula los valores de Shapley para los jugadores.
    
    Args:
        jugadores (list): Lista de jugadores (ej: ["A","B","C"]).
        valores_coaliciones (dict): Diccionario con los valores de las coaliciones.
                                 La clave es un frozenset de jugadores, ej:
                                 {frozenset(["A"]): 10, frozenset(["B","C"]): 25, ...}

    Devuelve:
        shapley (dict): Valor de Shapley de cada jugador.
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
                peso = math.factorial(len(S)) * math.factorial(n - len(S) - 1) / math.factorial(n)

                shapley[i] += peso * (v_S_i - v_S)

    return shapley