"""
Frontera sol-aire del C legacy (port para el harness de prueba de la Fase 4).

Reproduce literalmente `sol.h` + `Ta_Tc_DtaT` + `time_evolution_Tsa` del C
(`legacy_eh/2dTfree/`) para inyectar la MISMA ``Tsa(t)`` que el golden master y
poder validar el lazo temporal con tolerancia estricta. (En producción —Fase 5—
la frontera viene de ``Location.Tsa()`` con EPW+pvlib; esto es solo para validar.)
"""

import datetime
import math


# --- sol.h ---------------------------------------------------------------------

def _dia_juliano(dia, anio, mes):
    """juliano[m] = día del año (1-based) del `dia` de cada mes m=1..mes."""
    jul = [0] * (mes + 1)
    for m in range(1, mes + 1):
        jul[m] = datetime.date(anio, m, dia).timetuple().tm_yday
    return jul


def _solar_month(month, Lat, dia=15, anio=2011):
    """Devuelve (delta, D_dia, Ho) del mes, como en el `main` del C."""
    rad = math.pi / 180.0
    grad = 180.0 / math.pi
    jul = _dia_juliano(dia, anio, 12)
    delta = 23.45 * math.sin(rad * ((360.0 / 365.0) * (284.0 + jul[month])))
    orto = grad * math.acos((-1.0) * math.tan(rad * Lat) * math.tan(rad * delta))
    D_dia = (2.0 * orto) / 15.0          # duración del día (calculo_duracion_dia)
    Ho = 12.0 - (D_dia / 2.0)            # hora del amanecer (calculo_hora_orto)
    return delta, D_dia, Ho


# --- Ta_Tc_DtaT ----------------------------------------------------------------

def _y_factor(th, Ho, Hi):
    """Factor sinusoidal y(t) por tramos (amanecer/máximo), hora th en [0,24)."""
    pi = math.pi
    if th <= Ho:
        return (math.cos(pi * (Ho - th) / (24.0 + Ho - Hi)) + 1.0) / 2.0
    if th <= Hi:
        return (math.cos(pi * (th - Ho) / (Hi - Ho)) + 1.0) / 2.0
    return (math.cos(pi * (24.0 + Ho - th) / (24.0 + Ho - Hi)) + 1.0) / 2.0


def ta_tc_dtat(Tmax, Tmin, Hi, Ho):
    """Port de `Ta_Tc_DtaT`: temperatura media Ta, neutra Tc y semi-amplitud DtaT."""
    Ta = 0.0
    for t in range(0, 86401):          # t = 0..86400 inclusive (86401 muestras)
        y = _y_factor(t / 3600.0, Ho, Hi)
        Ta += y * Tmin + (1.0 - y) * Tmax
    Ta /= 86400.0                       # /86400 (no 86401) — fiel al C
    Tc = 0.54 * Ta + 13.5
    Delta = Tmax - Tmin
    bins = [(13, 2.5), (16, 3.0), (19, 3.5), (24, 4.0), (28, 4.5),
            (33, 5.0), (38, 5.5), (45, 6.0), (52, 6.5)]
    tmp2 = 7.0 / 2.0
    for hi_lim, val in bins:
        if Delta < hi_lim:
            tmp2 = val / 2.0
            break
    return Ta, Tc, tmp2


# --- time_evolution_Tsa --------------------------------------------------------

class Boundary:
    """Genera ``Tsa(t)``, ``Ta(t)`` e ``Is(t)`` idénticas al C para un mes/clima."""

    def __init__(self, params):
        g = lambda k: float(params[k])
        self.Ig, self.Id, self.Ib = g("Ig"), g("Id"), g("Ib")
        self.A = g("a")
        self.ho = g("ho")
        self.Tmax, self.Tmin = g("Tmax"), g("Tmin")
        self.Hi = g("tT")                 # hora de Tmax (t_Tamax)
        self.beta, self.gamma = g("beta"), g("gamma")
        self.Lat = g("Lat")
        month = int(params["mes"])
        delta, D_dia, Ho = _solar_month(month, self.Lat)
        self.delta = delta
        self.Ho = Ho
        self.tau = D_dia * 2.0            # el C pasa tau[month]*2
        self.Ta_mean, self.Tc, self.DtaT = ta_tc_dtat(
            self.Tmax, self.Tmin, self.Hi, self.Ho)

    def tsa(self, t, a=None):
        """Tsa sol-aire en el instante t [s]. Devuelve (Tsa, Ta_inst, Is)."""
        if a is None:
            a = self.A
        pi = math.pi
        X = pi / 180.0
        th = t / 3600.0
        y = _y_factor(th, self.Ho, self.Hi)
        tm2 = y * self.Tmin + (1.0 - y) * self.Tmax
        Ta_inst = tm2

        CF = 3.9 * (1.0 - self.beta / 90.0)
        if self.beta > 90.0:
            CF = 0.0

        Ibtheta = Idtheta = 0.0
        s = math.sin(2.0 * math.pi * (th / self.tau - self.Ho / self.tau))
        if s > 0.0:
            omega = th * 15.0 - 180.0
            phi, delta, beta, gamma = self.Lat, self.delta, self.beta, self.gamma
            thetaz = (math.cos(phi * X) * math.cos(delta * X) * math.cos(omega * X)
                      + math.sin(phi * X) * math.sin(delta * X))
            theta = (math.sin(delta * X) * math.sin(phi * X) * math.cos(beta * X)
                     - math.sin(delta * X) * math.cos(phi * X) * math.sin(beta * X) * math.cos(gamma * X)
                     + math.cos(delta * X) * math.cos(phi * X) * math.cos(beta * X) * math.cos(omega * X)
                     + math.cos(delta * X) * math.sin(phi * X) * math.sin(beta * X) * math.cos(gamma * X) * math.cos(omega * X)
                     + math.cos(delta * X) * math.sin(beta * X) * math.sin(gamma * X) * math.sin(omega * X))
            theta = math.acos(theta)
            thetaz = math.acos(thetaz)
            Ibtheta = self.Ib * s / math.cos(thetaz) * math.cos(theta)
            Idtheta = self.Id * s * (1.0 - beta / 180.0)
            if Ibtheta < 0.0:
                Ibtheta = 0.0

        Is = Idtheta + Ibtheta
        tsa = tm2 + a * (Idtheta + Ibtheta) / self.ho - CF
        return tsa, Ta_inst, Is
