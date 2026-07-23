"""
P2 — Factor de onda larga RF lineal en la inclinación.

Con ``absortance = 0`` la temperatura sol–aire se reduce a ``Tsa = Ta − RF``,
así que la diferencia ``Ta − Tsa`` mide RF directamente:

  RF(β) = 3.9·(1 − β/90°) °C  →  3.9 (techo), 1.95 (45°), 0 (muro), 0 (>90°).

En 0° y 90° es idéntico a la regla binaria previa (sin cambios en los casos
existentes); en inclinaciones intermedias elimina la discontinuidad y
restaura la regla lineal de la herramienta de 2016.

Correr con pytest o como script:
    .venv/bin/python tests/test_rf.py
"""

import os

import numpy as np

import enerhabitat as eh
from enerhabitat.config import config

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")


def _rf_at(tilt):
    config.file = os.path.join(HERE, "materials.ini")
    loc = eh.Location(EPW)
    loc.meanDay(month=5, year=2025)
    s = eh.System(loc, tilt=tilt, azimuth=0, absortance=0.0, layers=[("EPS", 0.1)])
    df = s.Tsa()
    rf = (df["Ta"] - df["Tsa"]).to_numpy()
    assert np.allclose(rf, rf[0], atol=1e-9), "RF debe ser constante en el día"
    return float(rf[0])


def test_rf_linear_in_tilt():
    expected = {0: 3.9, 45: 1.95, 90: 0.0, 120: 0.0}
    for tilt, rf in expected.items():
        got = _rf_at(tilt)
        assert abs(got - rf) <= 1e-9, f"tilt={tilt}: RF={got} != {rf}"


if __name__ == "__main__":
    test_rf_linear_in_tilt()
    print("PASS  test_rf_linear_in_tilt")
    print("\nP2: RF lineal 3.9→0 entre 0° y 90° ✅")
