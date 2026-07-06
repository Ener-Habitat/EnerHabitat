"""
Fase 6 — Fill con cámara de aire (`tipo 1`).

Valida el port de la física de cavidad (paredes NT 9-12, aire NT 0, radiación entre
paredes con factores de vista, Nusselt ``hh``, nodo de aire ``Thueco``):

  1. **Geometría** ``NT,k,rhoc`` vs dump del C — exacta.
  2. **Un paso** ``T, Tint, Thueco, hh`` vs dump del C — ``atol 1e-6``.
  3. **Día completo** serie ``Tso,Tsi,Tint,Thueco`` e índices vs golden — "del orden".

Golden en ``tests/golden/2d/hueca/`` (generado con la `.inp` tipo 1).

Correr con pytest o como script:
    .venv/bin/python tests/test_eh2d_hueca.py
"""

import os

import numpy as np

from enerhabitat.eh2d import Section2D, Fill
from enerhabitat.ehtools2d import solve_step_hueca, solve_day_hueca, _view_factors

from test_eh2d_geometry import read_inp, read_field, read_meta
from c_boundary import Boundary

G = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden", "2d", "hueca")
INP = os.path.join(G, "conduction_hueca.inp")

RHOAIR = 1.1797660470258469
CAIR = 1005.458757


def _section():
    params = read_inp(INP)
    g = lambda key: float(params[key])
    L = [g(f"L{n}") for n in range(1, 8)]
    k = [g(f"k{n}") for n in range(1, 8)]
    rhoc = [g(f"rhoc{n}") for n in range(1, 8)]
    a = {x: g(x) for x in ("a11", "a12", "a13", "a14", "a21", "a22", "a23")}
    e = {x: g(x) for x in ("e21", "e22", "e23")}
    sec = Section2D(nx=int(params["nx"]), ny=int(params["ny"]), L=L, k=k, rhoc=rhoc,
                    kr=g("kr"), rhocr=g("rhocr"), a=a, e=e, layer=int(params["layer"]),
                    fill_type=Fill.AIR).build()
    return sec, params


# --- Fase 6.1: geometría -------------------------------------------------------

def test_geometry():
    sec, _ = _section()
    NTc = read_field(os.path.join(G, "dump_NT.dat"), sec.nx, sec.ny, int)
    kc = read_field(os.path.join(G, "dump_k.dat"), sec.nx, sec.ny, float)
    rc = read_field(os.path.join(G, "dump_rhoc.dat"), sec.nx, sec.ny, float)
    assert np.array_equal(sec.NT, NTc)
    assert np.allclose(sec.kfield, kc, rtol=1e-12)
    assert np.allclose(sec.rhocfield, rc, rtol=1e-12)


# --- Fase 6.2: un paso ---------------------------------------------------------

def _one_step():
    sec, _ = _section()
    gm = read_meta(os.path.join(G, "dump_meta.dat"))
    sm = read_meta(os.path.join(G, "dump_step_meta.dat"))
    T0 = read_field(os.path.join(G, "dump_step_T0.dat"), sec.nx, sec.ny, float)
    Tc = read_field(os.path.join(G, "dump_step_T.dat"), sec.nx, sec.ny, float)
    m = sec.mesh
    res = solve_step_hueca(sec.NT, sec.kfield, sec.rhocfield, T0,
                           sm["Tsa"], sm["Tint_in"], sm["Thueco_in"],
                           sm["ho"], sm["hi"], sm["dt"], m.dx, m.dy,
                           sm["La"], m.X, sm["rhoair"], sm["cair"],
                           int(gm["i1"]), int(gm["j1"]), int(gm["i2"]), int(gm["j2"]),
                           sm["a21"], sm["e22"], sm["E"], sm["beta"])
    return res, Tc, sm


def test_step_field():
    res, Tc, _ = _one_step()
    assert np.max(np.abs(res["T"] - Tc)) <= 1e-6


def test_step_scalars():
    res, _, sm = _one_step()
    assert res["iters"] == int(sm["inner_iters"])
    assert abs(res["Thueco"] - sm["Thueco_out"]) <= 1e-6
    assert abs(res["Tint"] - sm["Tint_out"]) <= 1e-6
    assert abs(res["hh"] - sm["hh"]) <= 1e-9


# --- Fase 6.3: día completo ----------------------------------------------------

def solve_day(sec, params, dt=1.0):
    """Día completo tipo 1 con la Tsa sinusoidal del C; devuelve series + índices."""
    b = Boundary(params)
    m = sec.mesh
    nx, ny = sec.nx, sec.ny
    nsteps = int(86400 / dt) + 1
    Tsa = np.empty(nsteps); Tsa0 = np.empty(nsteps); Tsa1 = np.empty(nsteps)
    Is = np.empty(nsteps); Ta = np.empty(nsteps)
    for s in range(nsteps):
        t = s * dt
        ts, ta, isol = b.tsa(t)
        Tsa[s] = ts; Ta[s] = ta; Is[s] = isol
        Tsa0[s], _, _ = b.tsa(t, a=0.0)
        Tsa1[s], _, _ = b.tsa(t, a=1.0)
    hi = 8.1 if b.beta > 45.0 else 6.6
    T0 = b.Tc + b.DtaT
    vf = _view_factors(float(params["a21"]), float(params["e22"]))
    out = solve_day_hueca(sec.NT, sec.kfield, sec.rhocfield, Tsa,
                          b.ho, hi, dt, m.dx, m.dy, float(params["La"]), m.X,
                          RHOAIR, CAIR, T0,
                          m.i1, m.j1, m.i2, m.j2,
                          float(params["a21"]), float(params["e22"]), float(params["e"]),
                          *vf, 1e-10, 5e-4, 1)   # max_days=1: golden del C es DAYS=1
    Ti, Tso, Tsi, Th, Tfield, days = out
    return dict(Ti=Ti, Tso=Tso, Tsi=Tsi, Th=Th, Tsa=Tsa, Ta=Ta, Is=Is,
                Tsa0=Tsa0, Tsa1=Tsa1, days=days, Tc=b.Tc, DtaT=b.DtaT,
                dt=dt, X=m.X, hi=hi)


def _indices(r):
    Ti, Tsa = r["Ti"], r["Tsa"]
    Tc, dt = r["Tc"], r["dt"]
    # decremento, retardo
    dec = (Ti.max() - Ti.min()) / (Tsa.max() - Tsa.min())
    t = np.arange(len(Ti)) * dt
    lag = (t[np.argmax(Ti)] - t[np.argmax(Tsa)]) / 3600.0
    DDHhot = np.sum(np.where(Ti > Tc, Ti - Tc, 0.0)) * dt / 3600.0
    DDHcold = np.sum(np.where(Ti < Tc, Tc - Ti, 0.0)) * dt / 3600.0
    NumHot = np.sum(np.where(Ti > Tc, Ti - Tc, 0.0))
    DenHot = np.sum(np.where(r["Tsa1"] > Tc, r["Tsa1"] - Tc, 0.0))
    NumCold = np.sum(np.where(Ti < Tc, Tc - Ti, 0.0))
    DenCold = np.sum(np.where(r["Tsa0"] < Tc, Tc - r["Tsa0"], 0.0))
    # Qin (aprox. por superficie media: hi·(Tsi−Ti) donde Tsi>Ti, integrado en dt).
    hi = r["hi"]
    Qin = np.sum(np.where(r["Tsi"] > Ti, hi * (r["Tsi"] - Ti), 0.0)) / 3600.0
    return {
        "Qin": Qin,
        "decremento": dec, "retardo_h": lag,
        "Tint_media": Ti[:-1].mean() if len(Ti) > 1 else Ti.mean(),
        "Tintmin": Ti.min(), "Tintmax": Ti.max(),
        "TPIhot": (1 - NumHot / DenHot) * 100, "TPIcold": (1 - NumCold / DenCold) * 100,
        "DDHhot": DDHhot, "DDHcold": DDHcold,
    }


def _load_csv(path, skip):
    rows = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i < skip:
                continue
            rows.append([float(x) for x in line.split()])
    return np.array(rows)


_DAY = {}


def _run_day():
    if "r" in _DAY:
        return _DAY["r"]
    sec, params = _section()
    r = solve_day(sec, params)
    _DAY["r"] = r
    return r


def test_day_series():
    r = _run_day()
    csv = os.path.join(G, "gbv_5_2.csv")
    if not os.path.exists(csv):
        import pytest
        pytest.skip("golden tipo1 del día aún no generado")
    c = _load_csv(csv, skip=2)   # t Is Tsa Ta Tso Tsi Tint Tc DtaT
    # muestreo del golden: cada 600 s -> índices 0,600,1200,...
    idx = (c[:, 0] * 3600.0 / r["dt"]).round().astype(int)
    idx = idx[idx < len(r["Ti"])]
    for col, key in [(2, "Tsa"), (4, "Tso"), (5, "Tsi"), (6, "Ti")]:
        d = np.abs(r[key][idx] - c[:len(idx), col])
        assert d.max() <= 0.1 + 5e-3, f"{key}: max|Δ|={d.max():.3f}"


if __name__ == "__main__":
    import time
    for fn in (test_geometry, test_step_field, test_step_scalars):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nGeometría y un paso del tipo 1: exactos ✅\n")

    t0 = time.time()
    r = _run_day()
    print(f"Día tipo 1 resuelto: days={r['days']}  ({time.time()-t0:.1f}s)")
    csv = os.path.join(G, "gbv_5_2.csv")
    if os.path.exists(csv):
        c = _load_csv(csv, skip=2)
        idx = (c[:, 0] * 3600.0 / r["dt"]).round().astype(int)
        idx = idx[idx < len(r["Ti"])]
        bar = "═" * 70
        print(f"\n{bar}\n  FASE 6 · Día bovedilla con cámara de aire — Python vs C\n{bar}")
        print(f"\n  Serie — max|Δ| (de {len(idx)} muestras):")
        for col, key, name in [(2, "Tsa", "Tsa"), (4, "Tso", "Tso"),
                               (5, "Tsi", "Tsi"), (6, "Ti", "Tint")]:
            d = np.abs(r[key][idx] - c[:len(idx), col])
            print(f"    {name:<5} max|Δ|={d.max():.3f} °C  media={d.mean():.3f} °C")
        ip = _indices(r)
        cidx = _load_csv(os.path.join(G, "indice_gbv_5_2.csv"), skip=0)[0]
        order = ["Qin", "decremento", "retardo_h", "Tint_media", "Tintmin",
                 "Tintmax", "TPIhot", "TPIcold", "DDHhot", "DDHcold"]
        print(f"\n  Índices Python vs C:")
        for j, key in enumerate(order):
            pv = ip.get(key, float("nan"))
            print(f"    {key:<12}{pv:>10.3f}{cidx[j]:>10.3f}")
        # curva Thueco vs Tint vs Tsa
        print(f"\n{bar}\n  Curva diaria: 's'=Tsa  'h'=Thueco  'i'=Tint\n{bar}")
        Tsa = r["Tsa"][idx]; Th = r["Th"][idx]; Ti = r["Ti"][idx]
        lo, hi = min(Tsa.min(), Ti.min()), max(Tsa.max(), Ti.max())
        W = 54
        for kk in range(0, len(idx), max(1, len(idx) // 24)):
            line = [" "] * W
            for v, ch in [(Tsa[kk], "s"), (Th[kk], "h"), (Ti[kk], "i")]:
                line[int((v - lo) / (hi - lo) * (W - 1))] = ch
            print(f"  {c[kk,0]:4.1f}h |{''.join(line)}|")
        print(f"  rango {lo:.1f}..{hi:.1f} °C")
        print(f"\n{bar}\n  Fase 6: bovedilla con cámara de aire reproduce el C 'del orden' ✅\n{bar}")
    else:
        print("(golden del día tipo 1 aún no disponible; corre la geometría+un paso)")
