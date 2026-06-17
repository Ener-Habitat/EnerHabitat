"""
Fase 8b — Vigueta y bovedilla en techos (API `System2D` + `Slab`).

Losa de techo (`tilt=0`) con tres materiales sólidos (topping, vigueta en L,
bovedilla) y N cavidades iguales — de aire (`Fill.AIR`, Nusselt de techo
Rayleigh + radiación por hueco) o de relleno (`Fill.SOLID`). Se verifica
con la MISMA metodología que el 1D / 8a (config.file → Location → meanDay → Tsa →
solve), NO contra el C bit-a-bit:

  - **metodología/flujo**: `solve()` devuelve `pandas.Series` alineada a `Tsa()`;
  - **periodicidad**: converge día-a-día antes de `max_days`;
  - **balance de energía**: en régimen `Qin ≈ Qout`;
  - **física**: el factor de decremento del techo con AIRE supera al de relleno
    aislante (EPS) — el hueco transfiere más que el relleno aislante;
  - **capas antes/después**: acabado (L1) + Slab + yeso resuelven (capas ajenas
    al elemento se apilan por fuera);
  - **orientación**: `Slab` exige `tilt=0` (muro falla);
  - **inspector**: la sección muestra las N cavidades, la vigueta en L y los 3
    materiales.

Malla 2D reducida (`config2d`) para que la prueba corra rápido.

    .venv/bin/python tests/test_eh2d_slab.py
"""

import os

import numpy as np
import pandas as pd

import enerhabitat as eh
from enerhabitat.config import config2d
from enerhabitat.eh2d import SlabSection

HERE = os.path.dirname(os.path.abspath(__file__))
EPW = os.path.join(HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw")
MATERIALS = os.path.join(HERE, "materials_2d.ini")

# Geometría del paper (Fig. 2b), mm→m: vigueta en L (web/foot), 3 cavidades.
# topping=L2+L3=100, topping_cap=L2=50 → el alma de la L sube L3+L4+L5+L6=150 (no la tapa L2).
GEOM = {"web": 0.025, "foot": 0.025, "shoulder": 0.050, "n_cavities": 3,
        "cavity_width": 0.103, "topping": 0.100, "topping_cap": 0.050,
        "cover_top": 0.030, "cavity": 0.040, "cover_bottom": 0.030}


def _setup():
    eh.config.file = MATERIALS
    config2d.nx, config2d.ny = 48, 64     # malla chica: prueba de API, no de precisión
    config2d.max_days = 30
    loc = eh.Location(EPW)
    loc.meanDay(month=5, year=2025)
    return loc


def _roof(fill_type=eh.Fill.AIR, fill_material=None):
    loc = _setup()
    slab = eh.Slab("Concreto", fill_type=fill_type, block_material="Bovedilla",
                   topping_material="Concreto", fill_material=fill_material,
                   emissivity=0.9, geometry=GEOM)
    roof = eh.System2D(location=loc)
    roof.tilt = 0
    roof.absortance = 0.3
    roof.layers = [("Aplanado", 0.003), slab, ("Yeso", 0.015)]
    return roof


_CACHE = {}


def _solved(key, **kw):
    if key not in _CACHE:
        r = _roof(**kw)
        ti = r.solve()
        _CACHE[key] = (r, ti)
    return _CACHE[key]


def _decrement(r):
    df = r.solve_dataframe
    Tsa = df["Tsa"].to_numpy(); Ti = df["Ti"].to_numpy()
    return (Ti.max() - Ti.min()) / (Tsa.max() - Tsa.min())


# --- pruebas -------------------------------------------------------------------

def test_methodology_returns_series():
    r, ti = _solved("aire")
    assert isinstance(ti, pd.Series)
    assert len(ti) == len(r.Tsa())
    assert ti.index.equals(r.Tsa().index)
    assert np.isfinite(ti.to_numpy()).all()


def test_periodicity():
    r, _ = _solved("aire")
    assert 0 < r.days < config2d.max_days


def test_energy_balance():
    r, _ = _solved("aire")
    qin, qout = r.energy_transfer, r.Qout
    rel = abs(qin - qout) / max(abs(qin), 1e-9)
    assert rel <= 0.02, f"Qin={qin:.1f} Qout={qout:.1f} ({rel:.1%})"


def test_air_vs_insulating_fill():
    ra, _ = _solved("aire")
    rf, _ = _solved("eps", fill_type=eh.Fill.SOLID, fill_material="EPS")
    da, dfll = _decrement(ra), _decrement(rf)
    assert da > dfll, f"decremento aire={da:.3f} debe superar relleno EPS={dfll:.3f}"


def test_orientation_guard():
    loc = _setup()
    slab = eh.Slab("Concreto", block_material="Bovedilla", emissivity=0.9, geometry=GEOM)
    wall = eh.System2D(location=loc)
    wall.tilt = 90                              # vigueta y bovedilla NO es para muro
    wall.layers = [("Aplanado", 0.003), slab, ("Yeso", 0.015)]
    try:
        wall.solve()
    except ValueError:
        pass
    else:
        raise AssertionError("se esperaba ValueError por tilt≠0")


def test_inspector_geometry():
    r = _roof()
    sec = r.section()
    assert isinstance(sec, SlabSection)
    assert len(sec.cav_i1) == GEOM["n_cavities"]
    nts = set(np.unique(sec.NT).tolist())
    assert {0, 9, 10, 11, 12}.issubset(nts)         # aire + 4 paredes por hueco
    # tres materiales sólidos distintos (topping, vigueta, bovedilla) en el campo k
    ks = sorted({round(float(v), 4) for v in np.unique(sec.kfield)})
    assert len(ks) >= 3


def test_parallel_matches_serial():
    # El paralelo (prange) es opción; debe reproducir al serial (Jacobi por filas).
    # Opt-in (EH_TEST_PARALLEL=1): compilar el kernel paralelo tarda minutos, así que
    # no se corre en la suite por default. Malla mínima para que sea rápido.
    if not os.environ.get("EH_TEST_PARALLEL"):
        print("SKIP test_parallel_matches_serial (define EH_TEST_PARALLEL=1)")
        return
    r = _roof()
    config2d.nx, config2d.ny, config2d.max_days = 24, 32, 3
    config2d.parallel = True
    tip = r.solve().to_numpy()
    config2d.parallel = False
    tis = r.solve().to_numpy()
    dmax = float(np.abs(tip - tis).max())
    assert dmax <= 1e-6, f"paralelo vs serial max|Δ|={dmax:.2e}"


def test_l_shape_cap_height():
    # vigueta de k distinto al topping para poder distinguir el alma de la tapa.
    g = {"web": 0.025, "foot": 0.025, "shoulder": 0.050, "n_cav": 3,
         "cavity_width": 0.103, "topping": 0.100, "topping_cap": 0.050,
         "cover_top": 0.030, "cavity": 0.040, "cover_bottom": 0.030}
    nx, ny = 48, 64
    L = [0.0, 0.200, 0, 0, 0, 0, 0]
    sec = SlabSection(nx=nx, ny=ny, L=L, k=[0, 1.4, 0, 0, 0, 0, 0],
                      rhoc=[0, 2e6, 0, 0, 0, 0, 0], layer=2, geom=g,
                      k_topping=1.4, rc_topping=2e6, k_rib=2.0, rc_rib=2.2e6,
                      k_block=0.5, rc_block=1.1e6, emissivity=0.9, beta=0.0,
                      hollow=True).build()
    jet, jcap, jeb = sec.info["jet"], sec.info["jcap"], sec.info["jeb"]
    col = sec.kfield[0]                       # columna del alma (i=0)
    # tapa L2 [jet, jcap): topping (1.4), no vigueta
    assert np.allclose(col[jet:jcap], 1.4), "la tapa de topping L2 no debe ser vigueta"
    # alma [jcap, jeb): vigueta (2.0)
    assert np.allclose(col[jcap:jeb], 2.0), "el alma de la L debe ser vigueta de jcap a la base"
    # altura de la L ≈ L3+L4+L5+L6 = 150 mm (topping - topping_cap + resto)
    h_mm = (jeb - jcap) * sec.mesh.dy * 1000.0
    assert abs(h_mm - 150.0) < 2.0 * sec.mesh.dy * 1000.0, f"altura L={h_mm:.1f} mm (esperado ~150)"


def _demo():
    ra, _ = _solved("aire")
    rf, _ = _solved("eps", fill_type=eh.Fill.SOLID, fill_material="EPS")
    bar = "═" * 70
    print(bar)
    print("  FASE 8b · Techo de vigueta y bovedilla (System2D + Slab)")
    print(bar)
    ra.info()
    df = ra.solve_dataframe
    Tsa = df["Tsa"].to_numpy(); Ti = df["Ti"].to_numpy(); Th = df["Thueco"].to_numpy()
    print(f"\n  Convergió en {ra.days} días (tope {config2d.max_days})")
    qin, qout = ra.energy_transfer, ra.Qout
    rel = abs(qin - qout) / max(abs(qin), 1e-9)
    print(f"  Balance: Qin={qin:.1f}  Qout={qout:.1f}  desbalance={rel:.2%}")
    print(f"  Factor de decremento — AIRE: {_decrement(ra):.3f}   "
          f"RELLENA(EPS): {_decrement(rf):.3f}")
    print(f"  Tint(aire) rango {Ti.min():.2f}..{Ti.max():.2f} °C  ·  "
          f"Thueco {np.nanmin(Th):.2f}..{np.nanmax(Th):.2f} °C")

    print(f"\n{bar}\n  Sección a escala (techo, exterior arriba):\n{bar}")
    r2 = _roof()
    r2.section_report()
    print()
    r2.preview(field="materials", backend="ascii")

    print(f"\n{bar}\n  Curva diaria: 's'=Tsa  'h'=Thueco  'i'=Tint\n{bar}")
    n = len(Ti)
    lo, hi = min(Tsa.min(), Ti.min()), max(Tsa.max(), Ti.max())
    W = 54
    hrs = np.arange(n) * eh.config.dt / 3600.0
    for r in range(0, n, max(1, n // 24)):
        line = [" "] * W
        for v, ch in [(Tsa[r], "s"), (Th[r], "h"), (Ti[r], "i")]:
            line[int((v - lo) / (hi - lo) * (W - 1))] = ch
        print(f"  {hrs[r]:4.1f}h |{''.join(line)}|")
    print(f"  rango {lo:.1f}..{hi:.1f} °C")
    print(f"\n{bar}\n  Fase 8b: vigueta y bovedilla en techo, end-to-end ✅\n{bar}")


if __name__ == "__main__":
    for fn in (test_methodology_returns_series, test_periodicity, test_energy_balance,
               test_air_vs_insulating_fill, test_orientation_guard,
               test_inspector_geometry, test_l_shape_cap_height,
               test_parallel_matches_serial):
        fn()
        print(f"PASS  {fn.__name__}")
    print()
    _demo()
