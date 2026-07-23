"""
Propiedades del aire de las correlaciones de cavidad (techos, Hollands).

La ν = 1.11e-5 m²/s hardcodeada del C correspondía a aire a ~240 K,
inconsistente con el resto del conjunto (~300 K). Ahora se calcula:

  μ(300 K) por Sutherland  →  1.846e-5 Pa·s (Incropera A.4),
  ν = μ/ρ_aire  y  α = k/(ρ·c)  con la densidad CONFIGURABLE
  → ν y α termodinámicamente consistentes entre sí.

Correr con pytest o como script:
    .venv/bin/python tests/test_air_properties.py
"""

from enerhabitat.ehtools2d import _MU_AIR, _K_AIR, _BETA_EXP, _T_AIR_REF

RHOAIR = 1.1797660470258469
CAIR = 1005.458757


def test_sutherland_matches_incropera_300K():
    assert _T_AIR_REF == 300.0
    assert abs(_MU_AIR - 1.846e-5) <= 0.01e-5, f"mu={_MU_AIR:.4e}"


def test_nu_consistent_with_config_density():
    nu = _MU_AIR / RHOAIR
    # Incropera tabula 1.589e-5 con rho=1.1614; con la rho del paquete ~1.56e-5
    assert 1.5e-5 <= nu <= 1.62e-5, f"nu={nu:.4e}"
    # y ya no el valor de ~240 K heredado del C:
    assert abs(nu - 1.11e-5) > 3e-6


def test_alpha_beta_anchored_at_300K():
    alpha = _K_AIR / RHOAIR / CAIR
    assert 2.0e-5 <= alpha <= 2.4e-5           # Incropera 300 K: 2.25e-5
    assert abs(_BETA_EXP - 1.0 / 300.0) < 1e-12


if __name__ == "__main__":
    for fn in (test_sutherland_matches_incropera_300K,
               test_nu_consistent_with_config_density,
               test_alpha_beta_anchored_at_300K):
        fn()
        print(f"PASS  {fn.__name__}")
    nu = _MU_AIR / RHOAIR
    print(f"\nν efectiva = {nu:.4e} m²/s (antes 1.11e-5) ✅")
