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
    Compute the ambient temperature and add a 'Ta' column to the DataFrame.

    Args:
        df (pd.DataFrame): DataFrame whose 'index' represents the times.
        Tmin (float): Minimum temperature.
        Tmax (float): Maximum temperature.
        Ho (float): Sunrise hour (in hours).
        Hi (float): Hour of maximum temperature (in hours).

    Returns:
        pd.DataFrame: DataFrame with a new Ta column holding the ambient temperature.
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
        tmp2 = 0  # Optional, to cover any unhandled case, although the ranges above are exhaustive

    return tmp2

def get_sunrise_sunset_times(df):
    """
    Function to compute Ho and Hi.
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

def set_construction(materials, layers):
    """
    Update the cs dictionary with the material properties and the L values given in the tuples.

    Args:
        materials (dict): Dictionary with the material properties.
        layers (list): List of tuples, where each tuple holds the material and its L value.

    Returns:
        dict: Updated cs dictionary.
    """
    cs = {}
    for i, (material, L) in enumerate(layers, start=1):
        layer = f"L{i}"
        cs[layer] = {
            "L": L,
            "material": materials[material]
        }
    return cs

def get_total_L(cs):
    L_total = sum([cs[L]["L"] for L in cs.keys()])
    return L_total

def set_k_rhoc(cs, nx):
    """
    Map the physical layers onto the uniform 1D mesh (interface-aware).

    Assignment by cumulative coordinates:

    - ``k_array[i]``: conductivity of the material containing the cell
      **centre** (reference/reporting value; the solver uses ``Gf``).
    - ``rhoc_array[i]``: thickness-weighted average of ρc over the cell, so the
      total thermal mass ``Σ ρc_j·L_j`` is conserved exactly.
    - ``Gf[f]``: face conductance per unit area (W/m²K) between the centres of
      cells ``f`` and ``f+1``, from the exact series resistance
      ``∫ dx'/k(x')`` across the span. Material interfaces may fall anywhere
      inside a cell, and layers thinner than ``Δx`` (metal sheets, membranes)
      contribute their true resistance and mass at any position and any
      ``nx``. For a single-material span it reduces to ``k/Δx``; for an
      interface exactly at the face it reduces to the harmonic mean
      ``2·k_L·k_R/(k_L+k_R)/Δx``.

    Args:
        cs (dict): Dictionary with the constructive-system configuration.
        nx (int): Number of discretisation elements.

    Returns:
        tuple: (k_array, rhoc_array, dx, Gf) with ``Gf`` of shape ``(nx-1,)``.
    """
    L_total = get_total_L(cs)
    dx = L_total / nx

    n_lay = len(cs)
    bounds = np.zeros(n_lay + 1)
    k_lay = np.zeros(n_lay)
    rhoc_lay = np.zeros(n_lay)
    for j, L in enumerate(cs.keys()):
        bounds[j + 1] = bounds[j] + cs[L]['L']
        k_lay[j] = cs[L]['material'].k
        rhoc_lay[j] = cs[L]['material'].rho * cs[L]['material'].c
    bounds[-1] = L_total   # remove float accumulation dust

    def overlaps(x0, x1):
        """(length, layer index) of the intersection of [x0,x1] with each layer."""
        out = []
        for j in range(n_lay):
            lo = x0 if x0 > bounds[j] else bounds[j]
            hi = x1 if x1 < bounds[j + 1] else bounds[j + 1]
            if hi > lo:
                out.append((hi - lo, j))
        return out

    k_array = np.zeros(nx)
    rhoc_array = np.zeros(nx)
    for i in range(nx):
        x0 = i * dx
        x1 = x0 + dx
        xc = 0.5 * (x0 + x1)
        j = int(np.searchsorted(bounds, xc, side='right')) - 1
        j = 0 if j < 0 else (n_lay - 1 if j > n_lay - 1 else j)
        k_array[i] = k_lay[j]
        acc = 0.0
        for seg, jj in overlaps(x0, x1):
            acc += seg * rhoc_lay[jj]
        rhoc_array[i] = acc / dx

    Gf = np.zeros(max(nx - 1, 0))
    for f in range(nx - 1):
        x0 = (f + 0.5) * dx
        x1 = x0 + dx
        R = 0.0
        for seg, jj in overlaps(x0, x1):
            R += seg / k_lay[jj]
        Gf[f] = 1.0 / R

    return k_array, rhoc_array, dx, Gf

def prepare_static_coefficients(k_array, rhoc_array, dx, dt, ho, hi,
                                interface_cond=None):
    """
    Precompute mass and conductive coefficients that remain constant throughout the simulation.

    Args:
        k_array (numpy.ndarray): Thermal conductivity per node.
        rhoc_array (numpy.ndarray): Density × specific-heat product per node.
        dx (float): Control-volume size.
        dt (float): Time step.
        ho (float): Outdoor convective coefficient.
        hi (float): Indoor convective coefficient.
        interface_cond (numpy.ndarray, optional): Face conductances per unit
            area (W/m²K), shape (nx-1,) — e.g. the interface-aware ``Gf`` from
            :func:`set_k_rhoc`. If ``None``, falls back to the harmonic mean of
            the neighbouring cell conductivities (exact only when the material
            interface lies on the face).

    Returns:
        tuple: (mass_coeff, a_static, b_static, c_static) where:
            - mass_coeff (numpy.ndarray): Thermal-capacitance coefficients per node.
            - a_static (numpy.ndarray): Main diagonal of the tridiagonal system.
            - b_static (numpy.ndarray): Upper diagonal of the tridiagonal system.
            - c_static (numpy.ndarray): Lower diagonal of the tridiagonal system.
    """
    nx = k_array.shape[0]
    mass_coeff = rhoc_array * (dx / dt)

    if nx <= 1:
        a_static = np.empty(1, dtype=np.float64)
        b_static = np.zeros(1, dtype=np.float64)
        c_static = np.zeros(1, dtype=np.float64)
        a_static[0] = mass_coeff[0] + ho + hi
        return mass_coeff, a_static, b_static, c_static

    if interface_cond is None:
        inv_dx = 1.0 / dx
        interface_cond = np.empty(nx - 1, dtype=np.float64)
        for i in range(nx - 1):
            k_left = k_array[i]
            k_right = k_array[i + 1]
            interface_cond[i] = (2.0 * k_left * k_right) / (k_left + k_right) * inv_dx
    else:
        interface_cond = np.asarray(interface_cond, dtype=np.float64)

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
    Update in-place the right-hand-side vector of the tridiagonal system.

    Parameters:
        mass_coeff (numpy.ndarray): Precomputed thermal-capacitance coefficients per node.
        T (numpy.ndarray): Current domain temperatures.
        To (float): Outdoor temperature.
        ho (float): Outdoor convective coefficient.
        Ti (float): Indoor temperature.
        hi (float): Indoor convective coefficient.
        d (numpy.ndarray): Destination array for the thermal source term.
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
    Solve the equation system with the TDMA method and update the temperatures for the next time step.

    Args:
        a (numpy.ndarray): Array of a coefficients.
        b (numpy.ndarray): Array of b coefficients.
        c (numpy.ndarray): Array of c coefficients.
        d (numpy.ndarray): Array of d coefficients.
        T (numpy.ndarray): Temperature array.
        nx (int): Number of discretisation elements.
        Tint (float): Indoor temperature.
        capacitance_factor (float): Precomputed lumped-capacitance factor for the indoor space.
        P (numpy.ndarray): Auxiliary array for the forward-sweep phase.
        Q (numpy.ndarray): Auxiliary array for the forward-sweep phase.
        Tn (numpy.ndarray): Auxiliary array for the back substitution.

    Returns:
        tuple: (T, Tint) with the updated wall temperatures and the indoor temperature.
    """

    # Initialise P and Q
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
    Solve the tridiagonal system with TDMA for the air-conditioned mode, where the
    indoor temperature Tint is held constant (setpoint) instead of evolving with the
    lumped-capacitance model.

    Args:
        a (numpy.ndarray): Main diagonal of the tridiagonal system.
        b (numpy.ndarray): Upper diagonal of the tridiagonal system.
        c (numpy.ndarray): Lower diagonal of the tridiagonal system.
        d (numpy.ndarray): Right-hand-side vector (thermal source term).
        T (numpy.ndarray): Temperature array, updated in-place.
        nx (int): Number of discretisation elements.
        Tint (float): Indoor temperature (setpoint), returned unchanged.
        P (numpy.ndarray): Auxiliary array for the forward sweep.
        Q (numpy.ndarray): Auxiliary array for the forward sweep.
        Tn (numpy.ndarray): Auxiliary array for the back substitution.

    Returns:
        tuple: (T, Tint) with the updated wall temperatures and the indoor temperature
        (setpoint) unchanged.
    """

    # Initialise P and Q
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
