import numpy as np
import itertools
import math

def calcular_valores_coalicion(m, posiciones_genes, normalizar=True, metodo='cv'):
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
        suma_procesada = np.log2(1+suma_procesada)  # log(x + 1)
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

def calcula_valores_coaliciones(adata, nombres, indices):
    combinaciones_todos_tamanos = []

    for r in range(1, len(nombres) + 1):
        combinaciones_todos_tamanos.extend(itertools.combinations(indices, r))

    #for combinacion in combinaciones_todos_tamanos:
    #   print(combinacion)

    #with open("/Users/jose/Desktop/GameTheory_project/combinaciones.csv", "w", newline="") as archivo_csv:
    #    escritor = csv.writer(archivo_csv)
    #    escritor.writerows(combinaciones_todos_tamanos)

    print("Tamaño combinaciones de todos los tamaños:")
    print(len(combinaciones_todos_tamanos))

    # ESTRATEGIA HÍBRIDA:
    # 1. Matriz RAW para calcular_CV_grupo_genes (que normaliza internamente)
    matriz_raw = np.array(adata.layers["counts"].toarray()).T

    # 2. Matriz normalizada para cálculos individuales de CV, varianza y media
    matriz_norm = np.array(adata.X.toarray()).T

    num_filas, num_columnas = matriz_raw.shape

    print("NOTA: Estrategia híbrida implementada:")
    print("      - matriz_raw: Para calcular_CV_grupo_genes (normaliza internamente)")
    print("      - matriz_norm: Para cálculos individuales de CV, varianza y media")

    # Para mantener compatibilidad, usamos matriz_raw como 'matriz' principal
    matriz = matriz_raw

    print(f"Número de filas (genes): {num_filas}")
    print(f"Número de columnas (células): {num_columnas}")

    ganancias = {}

    # MÉTODO DE CÁLCULO: Elegir entre 'cv' (coeficiente de variación) o 'varianza'
    metodo_calculo = 'varianza'
    print(f"Método de cálculo seleccionado: {metodo_calculo}")

    for combinacion in combinaciones_todos_tamanos:
        
        ganancias[combinacion] = calcular_valores_coalicion(matriz_raw, np.array(combinacion), normalizar=True, metodo=metodo_calculo)

    print("Ganancias de todas las coaliciones:")

    i=0
    gan = np.zeros(len(ganancias))

    for combinacion, valor in ganancias.items():
        gan[i] = valor
        print(f"Combinación: {combinacion}, Valor: {valor}")
        i += 1
        
    return matriz_norm, ganancias


from combinations import ganancias_cython

def calcula_valores_coaliciones_cython(adata, nombres, indices):
    combinaciones_todos_tamanos = []
    
    max_cols = len(nombres)

    for r in range(1, len(nombres) + 1):
        curr_combination = itertools.combinations(indices, r)
        combinaciones_todos_tamanos.extend(curr_combination)        

    print("Tamaño combinaciones de todos los tamaños:")
    print(len(combinaciones_todos_tamanos))
    print("Max cols: ", max_cols)

    # ESTRATEGIA HÍBRIDA:
    # 1. Matriz RAW para calcular_CV_grupo_genes (que normaliza internamente)
    matriz_raw = np.array(adata.layers["counts"].toarray()).T

    # 2. Matriz normalizada para cálculos individuales de CV, varianza y media
    matriz_norm = np.array(adata.X.toarray()).T

    num_filas, num_columnas = matriz_raw.shape

    # Para mantener compatibilidad, usamos matriz_raw como 'matriz' principal
    matriz = matriz_raw

    print(f"Número de filas (genes): {num_filas}")
    print(f"Número de columnas (células): {num_columnas}")  

    ganancias = ganancias_cython(combinaciones_todos_tamanos, matriz, max_cols)
    
    N = len(combinaciones_todos_tamanos)
    ganancias_ret = {}
    i = 0
    for combinacion in combinaciones_todos_tamanos:
        ganancias_ret[combinacion] = ganancias[i]
        i += 1
        
    
    return matriz_norm, ganancias_ret

