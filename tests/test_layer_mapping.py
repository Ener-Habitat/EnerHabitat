"""
P0-05 — Mapeo de capas 1D a la malla (interface-aware).

Valida el mapeo por coordenadas acumuladas + conductancias de cara por
resistencias en serie que sustituye al conteo `int(L/dx)+1` heredado:

  1. **Caso de la revisión**: dos capas de 50 mm con Nx=10 → 5+5 celdas
     (el mapeo viejo daba 6+4).
  2. **Conservación de masa térmica**: Σ ρc_i·dx == Σ ρc_j·L_j exacta.
  3. **Resistencia de caras exacta**: Σ 1/Gf == ∫ dx'/k(x') entre el primer y
     el último centro de celda (integrador analítico independiente).
  4. **Lámina subcelda (zinc)**: aporta exactamente L/k a la resistencia, en
     cualquier posición interior y a cualquier Nx.
  5. **Simetría ante permutación**: invertir el orden de capas espeja los
     arreglos y conserva R total y masa.
  6. **Capa homogénea = mapeo legado**: Gf == k/dx bit a bit (garantiza que
     las comparaciones 2D→1D de una capa no cambian).
  7. **Refinamiento a nivel solver**: Ti(Nx=200) ≈ Ti(Nx=800) en un sistema
     multicapa con lámina metálica.

Correr con pytest o como script:
    .venv/bin/python tests/test_layer_mapping.py
"""

import os
from types import SimpleNamespace

import numpy as np

from enerhabitat.ehtools import set_construction, set_k_rhoc

HERE = os.path.dirname(os.path.abspath(__file__))

MATS = {
    "Zinc":     SimpleNamespace(k=110.0, rho=7140, c=390),
    "EPS":      SimpleNamespace(k=0.035, rho=30,   c=1400),
    "Concreto": SimpleNamespace(k=1.8,   rho=2400, c=900),
    "Yeso":     SimpleNamespace(k=0.30,  rho=800,  c=1090),
}


def _map(layers, nx):
    cs = set_construction(MATS, layers)
    return set_k_rhoc(cs, nx)


def _analytic_R(layers, a, b):
    """∫_a^b dx'/k(x') con un integrador independiente del código probado."""
    R, x0 = 0.0, 0.0
    for name, L in layers:
        x1 = x0 + L
        lo, hi = max(a, x0), min(b, x1)
        if hi > lo:
            R += (hi - lo) / MATS[name].k
        x0 = x1
    return R


def test_review_case_five_five():
    layers = [("EPS", 0.05), ("Concreto", 0.05)]
    k, rhoc, dx, Gf = _map(layers, 10)
    assert np.all(k[:5] == MATS["EPS"].k)
    assert np.all(k[5:] == MATS["Concreto"].k)


def test_mass_conservation():
    layers = [("Yeso", 0.0125), ("Concreto", 0.117), ("EPS", 0.0431),
              ("Zinc", 0.0005)]
    for nx in (7, 50, 200):
        k, rhoc, dx, Gf = _map(layers, nx)
        # malla de nodos-en-superficie: los nodos extremos poseen medio volumen
        vol = np.full(nx, dx)
        vol[0] = vol[-1] = 0.5 * dx
        total = (rhoc * vol).sum()
        exact = sum(MATS[n].rho * MATS[n].c * L for n, L in layers)
        assert abs(total - exact) <= 1e-9 * exact, f"masa nx={nx}"


def test_face_resistance_exact():
    layers = [("EPS", 0.033), ("Concreto", 0.1007), ("Yeso", 0.0125)]
    for nx in (11, 200):
        k, rhoc, dx, Gf = _map(layers, nx)
        R_faces = np.sum(1.0 / Gf)
        # los tramos nodo-a-nodo cubren [0, L] exacto: Σ 1/Gf = R total
        R_exact = _analytic_R(layers, 0.0, sum(L for _, L in layers))
        assert abs(R_faces - R_exact) <= 1e-12 * R_exact, f"R nx={nx}"


def test_subcell_zinc_resistance():
    """La lámina aporta exactamente L/k, en cualquier posición interior."""
    base = [("EPS", 0.05), ("Concreto", 0.10)]
    with_z1 = [("EPS", 0.05), ("Zinc", 0.0005), ("Concreto", 0.10)]
    with_z2 = [("EPS", 0.025), ("Zinc", 0.0005), ("EPS", 0.025),
               ("Concreto", 0.10)]
    R_zinc = 0.0005 / MATS["Zinc"].k
    for nx in (50, 200):
        for layers in (with_z1, with_z2):
            k, rhoc, dx, Gf = _map(layers, nx)
            R = np.sum(1.0 / Gf)
            R_exp = _analytic_R(layers, 0.0, sum(L for _, L in layers))
            assert abs(R - R_exp) <= 1e-12 * R_exp
            # y la contribución del zinc está incluida (no desaparece):
            layers_no = [(n, L) for n, L in layers if n != "Zinc"]
            # misma malla geométrica no aplica (L_total difiere); comparamos
            # la resistencia analítica total contra total sin zinc:
            R_tot = _analytic_R(layers, 0.0, sum(L for _, L in layers))
            R_tot_no = _analytic_R(layers_no, 0.0,
                                   sum(L for _, L in layers_no))
            assert abs((R_tot - R_tot_no) - R_zinc) <= 1e-15


def test_permutation_symmetric():
    fwd = [("EPS", 0.05), ("Concreto", 0.05)]
    rev = list(reversed(fwd))
    kf, rcf, dxf, Gff = _map(fwd, 10)
    kr, rcr, dxr, Gfr = _map(rev, 10)
    assert np.array_equal(kf, kr[::-1])
    assert np.allclose(rcf, rcr[::-1], rtol=1e-14)
    assert np.allclose(Gff, Gfr[::-1], rtol=1e-14)
    assert abs(np.sum(1 / Gff) - np.sum(1 / Gfr)) <= 1e-12


def test_single_layer_matches_legacy():
    k, rhoc, dx, Gf = _map([("EPS", 0.1)], 200)
    assert np.all(k == MATS["EPS"].k)
    assert np.allclose(rhoc, MATS["EPS"].rho * MATS["EPS"].c, rtol=1e-12)
    assert np.allclose(Gf, MATS["EPS"].k / dx, rtol=1e-12)


def test_solver_refinement():
    """Sistema multicapa con lámina: Ti converge con la malla."""
    import enerhabitat as eh
    from enerhabitat.config import config
    config.file = os.path.join(HERE, "materials_2d.ini")
    # añadir la lámina de zinc al vuelo
    config.materials["Zinc"] = SimpleNamespace(k=110.0, rho=7140, c=390)
    loc = eh.Location(os.path.join(
        HERE, "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw"))
    loc.meanDay(month=5, year=2025)
    layers = [("Zinc", 0.0005), ("EPS", 0.05), ("Concreto", 0.10)]

    def solve_with(nx):
        config.Nx = nx
        s = eh.System(loc, tilt=0, azimuth=0, absortance=0.8,
                      layers=list(layers))
        Ti = s.solve().to_numpy()
        assert s.converged
        return Ti

    try:
        Ti_a = solve_with(200)
        Ti_b = solve_with(800)
    finally:
        config.Nx = 200
    d = np.max(np.abs(Ti_a - Ti_b))
    assert d <= 0.05, f"refinamiento: max|ΔTi|={d:.4f} °C"


if __name__ == "__main__":
    for fn in (test_review_case_five_five, test_mass_conservation,
               test_face_resistance_exact, test_subcell_zinc_resistance,
               test_permutation_symmetric, test_single_layer_matches_legacy,
               test_solver_refinement):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nP0-05: mapeo interface-aware — masa y resistencia exactas ✅")
