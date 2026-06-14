import pandas as pd
import numpy as np
import math
from numba import njit
from dateutil.parser import parse

"""
=============================
        meanDay tools
=============================
"""

def add_temperature_model(df, Tmin, Tmax, Ho, Hi):
    """
    Calcula la temperatura ambiente y agrega una columna 'Ta' al DataFrame.

    Args:
        df (pd.DataFrame): DataFrame con la columna 'index' que representa los tiempos.
        Tmin (float): Temperatura mínima.
        Tmax (float): Temperatura máxima.
        Ho (float): Hora de amanecer (en horas).
        Hi (float): Hora de máxima temperatura (en horas).

    Returns:
        pd.DataFrame: DataFrame con una nueva columna Ta que contiene la temperatura ambiente.
    """
    Ho_sec = Ho * 3600
    Hi_sec = Hi * 3600
    day_hours = 24 * 3600
    times = pd.to_datetime(df.index)
    y = np.zeros(len(times))
    
    for i, t in enumerate(times):
        t_sec = t.hour * 3600 + t.minute * 60 + t.second
        if t_sec <= Ho_sec:
            y[i] = (math.cos(math.pi * (Ho_sec - t_sec) / (day_hours + Ho_sec - Hi_sec)) + 1) / 2
        elif Ho_sec < t_sec <= Hi_sec:
            y[i] = (math.cos(math.pi * (t_sec - Ho_sec) / (Hi_sec - Ho_sec)) + 1) / 2
        else:
            y[i] = (math.cos(math.pi * (day_hours + Ho_sec - t_sec) / (day_hours + Ho_sec - Hi_sec)) + 1) / 2

    Ta = Tmin + (Tmax - Tmin) * (1 - y)
    df['Ta'] = Ta
    return df

def calculate_tTmaxTminTmax(mes, epw):
    epw_mes = epw.loc[epw.index.month==int(mes)]
    hora_minutos = epw_mes.resample('D').To.idxmax()
    hora = hora_minutos.dt.hour
    minuto = hora_minutos.dt.minute
    tTmax = hora.mean() +  minuto.mean()/60
    Tmin =  epw_mes.resample('D').To.min().resample('ME').mean().iloc[0]
    Tmax =  epw_mes.resample('D').To.max().resample('ME').mean().iloc[0]
    
    return tTmax,Tmin,Tmax

def add_IgIbId_Tn(df, epw, mes, f1, f2, timezone):
    epw_mes = epw.loc[epw.index.month==int(mes)]
    Irr = epw_mes.groupby(by=epw_mes.index.hour)[['Ig','Id','Ib']].mean()
    tiempo = pd.date_range(start=f1, end=parse(f2), freq='1h',tz=timezone)
    Irr.index = tiempo
    Irr = Irr.resample('1s').interpolate(method='time')
    df['Ig'] = Irr.Ig
    df['Ib'] = Irr.Ib
    df['Id'] = Irr.Id
    df.ffill(inplace=True)
    df['Tn'] = 13.5 + 0.54*df.Ta.mean()
    
    return df

@njit
def calculate_DtaTn(Delta):
    if Delta < 13:
        tmp2 = 2.5 / 2
    elif 13 <= Delta < 16:
        tmp2 = 3.0 / 2
    elif 16 <= Delta < 19:
        tmp2 = 3.5 / 2
    elif 19 <= Delta < 24:
        tmp2 = 4.0 / 2
    elif 24 <= Delta < 28:
        tmp2 = 4.5 / 2
    elif 28 <= Delta < 33:
        tmp2 = 5.0 / 2
    elif 33 <= Delta < 38:
        tmp2 = 5.5 / 2
    elif 38 <= Delta < 45:
        tmp2 = 6.0 / 2
    elif 45 <= Delta < 52:
        tmp2 = 6.5 / 2
    elif Delta >= 52:
        tmp2 = 7.0 / 2
    else:
        tmp2 = 0  # Opcional, para cubrir cualquier caso no contemplado, aunque el rango anterior es exhaustivo

    return tmp2

def get_sunrise_sunset_times(df):
    """
    Función para calcular Ho y Hi
    """
    sunrise_time = df[df['elevation'] >= 0].index[0]
    sunset_time = df[df['elevation'] >= 0].index[-1]
    
    Ho = sunrise_time.hour + sunrise_time.minute / 60
    Hi = sunset_time.hour + sunset_time.minute / 60
    
    return Ho, Hi

"""
=============================
        solveCS tools
=============================
"""

def set_construction(propiedades, tuplas):
    """
    Actualiza el diccionario cs con  las propiedades del material y los valores de L proporcionados en las tuplas.
    
    Argss:
        propiedades (dict): Diccionario con las propiedades de los materiales.
        tuplas (list): Lista de tuplas, donde cada tupla contiene el material y el valor de L.
    
    Returns:
        dict: Diccionario actualizado cs.
    """
    cs ={}
    for i, (material, L) in enumerate(tuplas, start=1):
        capa = f"L{i}"
        cs[capa] = {
            "L": L,
            "material": propiedades[material]
        }
    return cs

def get_total_L(cs):
    L_total = sum([cs[L]["L"] for L in cs.keys()])
    return L_total

def set_k_rhoc(cs, nx):
    """
    Calcula los arreglos de conductividad y el producto de calor específico y densidad
    para cada volumen de control, y también calcula el tamaño de cada volumen de control (dx).
    
    Args:
        cs (dict): Diccionario con la configuración del sistema constructivo.
        nx (int): Número de elementos de discretización.
    
    Returns:
        tuple : [ k_array, rhoc_array, dx ] donde k_array es el arreglo de conductividad,
        rhoc_array es el arreglo del producto de calor específico y densidad,
        y dx es el tamaño de cada volumen de control.
    """
    L_total = get_total_L(cs)
    dx = L_total / nx

    k_array = np.zeros(nx)
    rhoc_array = np.zeros(nx)

    # Inicializar la posición actual en el arreglo
    i = 0

    for L in cs.keys():
        L_value = cs[L]['L']
        k_value = cs[L]['material'].k
        rhoc_value = cs[L]['material'].rho * cs[L]['material'].c

        num_elements = int(L_value / dx)
        
        for j in range(num_elements):
            if i >= nx:
                break
            k_array[i] = k_value
            rhoc_array[i] = rhoc_value
            i += 1

        # Considerar promedio armónico solo con el primer vecino
        if i < nx and i > 0:
            k_array[i] = 2 * (k_array[i-1] * k_value) / (k_array[i-1] + k_value)
            rhoc_array[i] = rhoc_value
            i += 1

    return k_array, rhoc_array, dx

def prepare_static_coefficients(k_array, rhoc_array, dx, dt, ho, hi):
    """
    Precompute mass and conductive coefficients that remain constant throughout the simulation.

    Args:
        k_array (numpy.ndarray): Conductividad térmica por nodo.
        rhoc_array (numpy.ndarray): Producto densidad * calor específico por nodo.
        dx (float): Tamaño del volumen de control.
        dt (float): Paso de tiempo.
        ho (float): Coeficiente convectivo exterior.
        hi (float): Coeficiente convectivo interior.

    Returns:
        tuple: (mass_coeff, a_static, b_static, c_static) donde:
            - mass_coeff (numpy.ndarray): Coeficientes de capacidad térmica por nodo.
            - a_static (numpy.ndarray): Diagonal principal del sistema tridiagonal.
            - b_static (numpy.ndarray): Diagonal superior del sistema tridiagonal.
            - c_static (numpy.ndarray): Diagonal inferior del sistema tridiagonal.
    """
    nx = k_array.shape[0]
    mass_coeff = rhoc_array * (dx / dt)

    if nx <= 1:
        a_static = np.empty(1, dtype=np.float64)
        b_static = np.zeros(1, dtype=np.float64)
        c_static = np.zeros(1, dtype=np.float64)
        a_static[0] = mass_coeff[0] + ho + hi
        return mass_coeff, a_static, b_static, c_static

    inv_dx = 1.0 / dx
    interface_cond = np.empty(nx - 1, dtype=np.float64)

    for i in range(nx - 1):
        k_left = k_array[i]
        k_right = k_array[i + 1]
        interface_cond[i] = (2.0 * k_left * k_right) / (k_left + k_right) * inv_dx

    a_static = np.empty(nx, dtype=np.float64)
    b_static = np.empty(nx, dtype=np.float64)
    c_static = np.empty(nx, dtype=np.float64)

    cond_right = interface_cond[0]
    mass0 = mass_coeff[0]
    a_static[0] = mass0 + ho + cond_right
    b_static[0] = cond_right
    c_static[0] = 0.0

    for i in range(1, nx - 1):
        cond_left = interface_cond[i - 1]
        cond_right = interface_cond[i]
        mass_i = mass_coeff[i]

        a_static[i] = mass_i + cond_left + cond_right
        b_static[i] = cond_right
        c_static[i] = cond_left

    cond_left = interface_cond[nx - 2]
    mass_last = mass_coeff[nx - 1]
    a_static[nx - 1] = mass_last + cond_left + hi
    b_static[nx - 1] = 0.0
    c_static[nx - 1] = cond_left

    return mass_coeff, a_static, b_static, c_static

@njit(cache=True)
def calculate_coefficients(mass_coeff, T, To, ho, Ti, hi, d):
    """
    Actualiza in-place el vector de términos independientes del sistema tridiagonal.

    Parameters:
        mass_coeff (numpy.ndarray): Coeficientes de capacidad térmica precomputados por nodo.
        T (numpy.ndarray): Temperaturas actuales del dominio.
        To (float): Temperatura en el exterior.
        ho (float): Coeficiente convectivo exterior.
        Ti (float): Temperatura en el interior.
        hi (float): Coeficiente convectivo interior.
        d (numpy.ndarray): Arreglo destino para la fuente térmica.
    """
    nx = mass_coeff.shape[0]

    if nx == 1:
        mass0 = mass_coeff[0]
        d[0] = mass0 * T[0] + ho * To + hi * Ti
        return

    mass0 = mass_coeff[0]
    d[0] = mass0 * T[0] + ho * To

    for i in range(1, nx - 1):
        mass_i = mass_coeff[i]
        d[i] = mass_i * T[i]

    mass_last = mass_coeff[nx - 1]
    d[nx - 1] = mass_last * T[nx - 1] + hi * Ti

@njit(cache=True)
def solve_PQ(a, b, c, d, T, nx, Tint, capacitance_factor, P, Q, Tn):
    """
    Resuelve el sistema de ecuaciones usando el método TDMA y actualiza las temperaturas para el siguiente paso temporal.

    Args:
        a (numpy.ndarray): Arreglo de coeficientes a.
        b (numpy.ndarray): Arreglo de coeficientes b.
        c (numpy.ndarray): Arreglo de coeficientes c.
        d (numpy.ndarray): Arreglo de coeficientes d.
        T (numpy.ndarray): Arreglo de temperaturas.
        nx (int): Número de elementos de discretización.
        Tint (float): Temperatura interna.
        capacitance_factor (float): Factor lumped-capacitance precomputado para el recinto.
        P (numpy.ndarray): Arreglo auxiliar para la fase de forward sweep.
        Q (numpy.ndarray): Arreglo auxiliar para la fase de forward sweep.
        Tn (numpy.ndarray): Arreglo auxiliar para el back substitution.

    Returns:
        tuple: (T, Tint) con las temperaturas de muro actualizadas y la temperatura interior.
    """

    # Inicializar P y Q
    inv_a0 = 1.0 / a[0]
    P[0] = b[0] * inv_a0
    Q[0] = d[0] * inv_a0

    for i in range(1, nx):
        denom = a[i] - c[i] * P[i - 1]
        inv_denom = 1.0 / denom
        P[i] = b[i] * inv_denom
        Q[i] = (d[i] + c[i] * Q[i - 1]) * inv_denom

    Tn[nx - 1] = Q[nx - 1]
    for i in range(nx - 2, -1, -1):
        Tn[i] = P[i] * Tn[i + 1] + Q[i]

    for i in range(nx):
        T[i] = Tn[i]

    Tint = Tint + capacitance_factor * (T[nx - 1] - Tint)

    return T, Tint

@njit(cache=True)
def solve_PQ_AC(a, b, c, d, T, nx, Tint, P, Q, Tn):
    """
    Resuelve el sistema tridiagonal con TDMA para el modo con aire acondicionado,
    donde la temperatura interior Tint se mantiene constante (setpoint) en lugar
    de evolucionar con el modelo lumped-capacitance.

    Args:
        a (numpy.ndarray): Diagonal principal del sistema tridiagonal.
        b (numpy.ndarray): Diagonal superior del sistema tridiagonal.
        c (numpy.ndarray): Diagonal inferior del sistema tridiagonal.
        d (numpy.ndarray): Vector de términos independientes (fuente térmica).
        T (numpy.ndarray): Arreglo de temperaturas, actualizado in-place.
        nx (int): Número de elementos de discretización.
        Tint (float): Temperatura interior (setpoint), devuelta sin cambios.
        P (numpy.ndarray): Arreglo auxiliar para el forward sweep.
        Q (numpy.ndarray): Arreglo auxiliar para el forward sweep.
        Tn (numpy.ndarray): Arreglo auxiliar para el back substitution.

    Returns:
        tuple: (T, Tint) con las temperaturas de muro actualizadas y la temperatura
        interior (setpoint) sin cambios.
    """

    # Inicializar P y Q
    inv_a0 = 1.0 / a[0]
    P[0] = b[0] * inv_a0
    Q[0] = d[0] * inv_a0

    for i in range(1, nx):
        denom = a[i] - c[i] * P[i - 1]
        inv_denom = 1.0 / denom
        P[i] = b[i] * inv_denom
        Q[i] = (d[i] + c[i] * Q[i - 1]) * inv_denom

    Tn[nx - 1] = Q[nx - 1]
    for i in range(nx - 2, -1, -1):
        Tn[i] = P[i] * Tn[i + 1] + Q[i]

    for i in range(nx):
        T[i] = Tn[i]

    return T, Tint
