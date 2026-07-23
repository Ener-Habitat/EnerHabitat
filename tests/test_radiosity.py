"""
P0-02 — Radiosidad de cavidades (factores de transferencia de Gebhart).

Valida que los factores 𝔉 = ε²(I−(1−ε)F)⁻¹F usados en producción resuelven el
recinto gris difuso de 4 superficies exactamente:

  1. **Reciprocidad**: A_m·𝔉_mn == A_n·𝔉_nm.
  2. **Cierre**: Σ_n 𝔉_mn == ε sobre la matriz completa (incluida la diagonal
     𝔉_mm — radiación que vuelve a la propia superficie por reflexiones; en la
     forma pareada esa diagonal se cancela con T⁴−T⁴=0 y no se exporta).
  3. **Límite negro**: ε = 1 → 𝔉 == F (recupera el modelo del C).
  4. **Placas paralelas**: cavidad muy ancha → 𝔉_ud → 1/(2/ε − 1) analítico.
  5. **Conservación**: Σ_m A_m·q_m == 0 para temperaturas arbitrarias.
  6. **Dirección del cambio**: con ε<1 el intercambio del par dominante es
     MENOR que el pareado ε·F (el modelo viejo sobreestimaba).

Correr con pytest o como script:
    .venv/bin/python tests/test_radiosity.py
"""

import numpy as np

from enerhabitat.ehtools2d import (_view_factors, _view_factor_matrix,
                                   _transfer_factors)


def _G_matrix(a21, e22, eps):
    """Reconstruye la matriz 4×4 de 𝔉 desde la tupla de 12 nombres."""
    (Fud, Ful, Fur, Fru, Frd, Frl,
     Fdl, Fdr, Fdu, Flu, Flr, Fld) = _transfer_factors(a21, e22, eps)
    G = np.zeros((4, 4))   # orden u, d, l, r
    G[0, 1], G[0, 2], G[0, 3] = Fud, Ful, Fur
    G[1, 0], G[1, 2], G[1, 3] = Fdu, Fdl, Fdr
    G[2, 0], G[2, 3], G[2, 1] = Flu, Flr, Fld
    G[3, 0], G[3, 1], G[3, 2] = Fru, Frd, Frl
    return G


# geometría del bloque documentado: cavidad 0.16 × 0.08 m
A21, E22 = 0.16, 0.08
AREAS = np.array([A21, A21, E22, E22])   # u, d, l, r


def test_reciprocity():
    for eps in (0.3, 0.9):
        G = _G_matrix(A21, E22, eps)
        AG = AREAS[:, None] * G
        assert np.allclose(AG, AG.T, rtol=1e-12), f"reciprocidad ε={eps}"


def _G_full(a21, e22, eps):
    """𝔉 completa (con diagonal), calculada en el test desde F."""
    F = _view_factor_matrix(a21, e22)
    return (eps * eps) * np.linalg.solve(np.eye(4) - (1.0 - eps) * F, F)


def test_closure():
    for eps in (0.3, 0.9):
        Gfull = _G_full(A21, E22, eps)
        assert np.allclose(Gfull.sum(axis=1), eps, rtol=1e-12), f"cierre ε={eps}"
        # y la tupla exportada coincide con los fuera-de-diagonal:
        G = _G_matrix(A21, E22, eps)
        off = Gfull.copy(); np.fill_diagonal(off, 0.0)
        assert np.allclose(G, off, rtol=1e-12)


def test_black_limit_recovers_legacy():
    G = _G_matrix(A21, E22, 1.0)
    F = _view_factor_matrix(A21, E22)
    assert np.allclose(G, F, rtol=1e-12)


def test_parallel_plate_limit():
    eps = 0.9
    G = _G_matrix(1000.0, 1.0, eps)   # a21 >> e22: u y d dominan
    exact = 1.0 / (2.0 / eps - 1.0)   # = 0.8181...
    assert abs(G[0, 1] - exact) <= 1e-2 * exact, \
        f"placas paralelas: {G[0,1]:.5f} vs {exact:.5f}"


def test_energy_conservation():
    eps = 0.9
    G = _G_matrix(A21, E22, eps)
    T4 = (np.array([305.0, 295.0, 301.0, 288.0])) ** 4
    q = np.array([np.sum(G[m] * (T4[m] - T4)) for m in range(4)])
    total = np.sum(AREAS * q)
    scale = np.sum(AREAS * np.abs(q))
    assert abs(total) <= 1e-12 * scale, f"no conserva: {total:.3e}"


def test_grey_exchange_below_pairwise():
    eps = 0.9
    G = _G_matrix(A21, E22, eps)
    F = _view_factor_matrix(A21, E22)
    assert G[0, 1] < eps * F[0, 1], \
        "el par dominante debe intercambiar menos que el modelo pareado ε·F"


if __name__ == "__main__":
    for fn in (test_reciprocity, test_closure, test_black_limit_recovers_legacy,
               test_parallel_plate_limit, test_energy_conservation,
               test_grey_exchange_below_pairwise):
        fn()
        print(f"PASS  {fn.__name__}")
    print("\nP0-02: radiosidad (Gebhart) — reciprocidad, cierre y límites ✅")
