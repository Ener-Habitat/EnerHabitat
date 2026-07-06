import marimo

__generated_with = "0.23.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.md(r"""
    # EnerHabitat — demo / pruebas: techo de 10 cm de EPS

    Notebook de prueba del uso de la librería **EnerHabitat** sobre un sistema
    constructivo simple: una azotea aislada con **10 cm de EPS**
    (`k = 0.035 W/m·K`, `rho = 35 kg/m³`, `c = 1000 J/kg·K`),
    clima de **Cuernavaca, Morelos** (EPW TMYx), día promedio de **mayo**,
    superficie clara (`absortancia = 0.3`).

    No requiere instalar marimo en el proyecto: ábrelo con
    `uv run --with marimo marimo edit tests/eps_demo.py`.
    Las celdas con `assert` al final actúan como pruebas (smoke tests).
    """)
    return


@app.cell
def _():
    from pathlib import Path

    import numpy as np
    import pandas as pd

    import enerhabitat as eh

    HERE = Path(__file__).parent
    EPW = HERE / "MEX_MOR_Cuernavaca-Matamoros.Intl.AP.767260_TMYx.2004-2018.epw"
    MATERIALS = HERE / "materials.ini"

    MONTH = 5
    YEAR = 2023
    return EPW, MATERIALS, MONTH, YEAR, eh, np, pd


@app.cell
def _(MATERIALS, eh, mo):
    # 1) Materiales (requerido — el paquete no trae materiales por defecto)
    eh.config.file = str(MATERIALS)

    _eps = eh.config.materials["EPS"]
    mo.md(
        f"""
        ### Configuración global

        Materiales cargados: `{eh.config.materials_list()}`

        **EPS** → k = {_eps.k} W/m·K · rho = {_eps.rho} kg/m³ · c = {_eps.c} J/kg·K

        Discretización: Nx = {eh.config.Nx} volúmenes · dt = {eh.config.dt} s ·
        ho = {eh.config.ho} · hi = {eh.config.hi} W/m²·K
        """
    )
    return


@app.cell
def _(EPW, eh, mo):
    # 2) Ubicación a partir del EPW
    loc = eh.Location(str(EPW))
    mo.md(
        f"""
        ### Ubicación

        **{loc.city}** — lat {loc.latitude}° · lon {loc.longitude}° ·
        alt {loc.altitude} m · tz `{loc.timezone}`
        """
    )
    return (loc,)


@app.cell
def _(eh, loc):
    # 3) Sistema constructivo: azotea (tilt=0) de 10 cm de EPS, color claro
    techo = eh.System(loc)
    techo.tilt = 0          # 0 = techo horizontal (RF = 3.9 °C)
    techo.azimuth = 0
    techo.absortance = 0.3  # superficie clara
    techo.layers = [("EPS", 0.10)]
    return (techo,)


@app.cell
def _(MONTH, YEAR, loc, techo):
    # 4) Día promedio + temperatura sol-aire
    loc.meanDay(month=MONTH, year=YEAR)
    tsa = techo.Tsa()
    return (tsa,)


@app.cell
def _(mo, tsa):
    mo.md(f"""
    ### Día promedio (mayo) y temperatura sol-aire

    | Variable | mín | máx | media |
    | --- | --- | --- | --- |
    | Ta (ambiente)   | {tsa.Ta.min():.1f} | {tsa.Ta.max():.1f} | {tsa.Ta.mean():.1f} |
    | Tsa (sol-aire)  | {tsa.Tsa.min():.1f} | {tsa.Tsa.max():.1f} | {tsa.Tsa.mean():.1f} |
    | Is (irrad. plano) | {tsa.Is.min():.0f} | {tsa.Is.max():.0f} | {tsa.Is.mean():.0f} |
    | Tn (neutra)     | — | — | {tsa.Tn.mean():.1f} |

    De noche (`Is = 0`) la azotea pierde calor al cielo: `Tsa = Ta − 3.9 °C`.
    """)
    return


@app.cell
def _(pd, techo, tsa):
    # 5) Solución en flotación libre (sin aire acondicionado)
    ti = techo.solve()
    energy_transfer = techo.energy_transfer  # J/(m²·día)

    # DataFrame combinado para inspección/gráfica (Tsa muestreada a 10 min)
    serie = pd.concat([ti, tsa[["Ta", "Tsa", "Tn"]].asfreq("10min")], axis=1)
    return energy_transfer, serie, ti


@app.cell
def _(energy_transfer, mo, ti):
    mo.md(f"""
    ### Flotación libre — `solve()`

    Temperatura interior `Ti`: {ti.min():.1f} – {ti.max():.1f} °C
    (oscila mucho porque el EPS casi no tiene masa térmica).

    Energía transferida al interior: **{energy_transfer:.0f} J/(m²·día)**
    = {energy_transfer/3600:.2f} Wh/(m²·día).
    """)
    return


@app.cell
def _(serie):
    # Tabla de la serie combinada (marimo la renderiza como tabla interactiva)
    serie
    return


@app.cell
def _(mo, techo):
    # 6) Solución con aire acondicionado (setpoint de confort)
    techo.solveAC()
    cooling = techo.cooling_energy  # J/(m²·día)
    heating = techo.heating_energy  # J/(m²·día)

    mo.md(
        f"""
        ### Con aire acondicionado — `solveAC()`

        Demanda para mantener el setpoint de confort en un día promedio de mayo:

        | | J/(m²·día) | Wh/(m²·día) |
        | --- | --- | --- |
        | Enfriamiento | {cooling:,.0f} | {cooling/3600:.2f} |
        | Calefacción  | {heating:,.0f} | {heating/3600:.2f} |
        """
    )
    return cooling, heating


@app.cell
def _(mo):
    mo.md(r"""
    ## Pruebas (smoke tests)
    """)
    return


@app.cell
def _(cooling, eh, energy_transfer, heating, loc, np, ti, tsa):
    # --- Asserts: si alguno falla, marimo marca la celda con error ---
    checks = []

    def check(name, cond):
        checks.append((name, bool(cond)))
        assert cond, f"FALLÓ: {name}"

    # Materiales y ubicación
    check("EPS cargado con k=0.035", eh.config.materials["EPS"].k == 0.035)
    check("EPW leído (Cuernavaca, lat≈18.8)", "Cuernavaca" in loc.city and abs(loc.latitude - 18.835) < 1e-3)

    # Día promedio
    check("día promedio sin NaN", not tsa[["Ta", "Ig", "Ib", "Id", "Tsa", "Is"]].isna().any().any())
    check("serie de ~1 día (>80k filas a 1 s)", len(tsa) > 80_000)

    # Temperatura sol-aire
    night = tsa[tsa.Is == 0]
    check("RF de techo: Tsa = Ta - 3.9 de noche", np.allclose((night.Tsa - night.Ta), -3.9, atol=1e-6))
    check("ganancia solar: Tsa_max > Ta_max", tsa.Tsa.max() > tsa.Ta.max())

    # Flotación libre
    check("Ti sin NaN", not ti.isna().any())
    check("Ti dentro de cota física [Ta_min-4, Tsa_max]", ti.min() >= tsa.Ta.min() - 4 and ti.max() <= tsa.Tsa.max() + 1e-6)
    check("estado periódico: Ti[0] ≈ Ti[-1]", abs(ti.iloc[0] - ti.iloc[-1]) < 0.5)
    # energy_transfer = calor que ENTRA al recinto en flotación libre,
    # acumulado en el día promedio: Σ hi·(T_sup_int − Ti)·dt cuando T_sup_int > Ti
    check("energía transferida (la que entra) ≥ 0", energy_transfer >= 0)

    # Aire acondicionado
    check("enfriamiento ≥ 0", cooling >= 0)
    check("calefacción ≥ 0", heating >= 0)
    return (checks,)


@app.cell
def _(checks, mo):
    _rows = "\n".join(f"| {'✅' if ok else '❌'} | {name} |" for name, ok in checks)
    _npass = sum(ok for _, ok in checks)
    mo.md(
        f"""
        ### Resultado: {_npass}/{len(checks)} pruebas OK

        | | prueba |
        | --- | --- |
        {_rows}
        """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Prueba unitaria del kernel TDMA (no de los modos)

    `solve_PQ` (flotación libre) y `solve_PQ_AC` (A/C) comparten **el mismo**
    algoritmo TDMA, que solo resuelve el sistema tridiagonal del muro
    `[A]{T} = {d}`. Lo que distingue a los dos modos **no está en el TDMA**, sino
    en lo que pasa *después*: en flotación libre la temperatura interior evoluciona
    (`Ti += factor·(T_muro − Ti)`), mientras que con A/C se mantiene fija en el
    setpoint y el *caller* contabiliza la energía de enfriar/calentar.

    Por eso esta prueba **no** afirma que ambos modos den el mismo resultado físico
    (no lo dan: en una simulación real reciben un `d` distinto y divergen). Lo que
    hace es darles a las dos funciones **exactamente los mismos** coeficientes
    `a, b, c, d` y verificar que devuelven el mismo vector de temperaturas de muro.
    Es un *regression guard* del refactor que compiló `solve_PQ_AC` con numba:
    confirma que el kernel nuevo resuelve idéntico al kernel ya probado.
    (Se pasa `capacitance_factor = 0` a `solve_PQ` para anular su actualización de
    `Ti` y comparar únicamente el muro.)
    """)
    return


@app.cell
def _(np):
    from enerhabitat.ehtools import solve_PQ, solve_PQ_AC

    _nx = 8
    _a = np.full(_nx, 4.0)
    _b = np.array([1.0, 1, 1, 1, 1, 1, 1, 0])
    _c = np.array([0.0, 1, 1, 1, 1, 1, 1, 1])
    _d = np.array([10.0, 5, 5, 5, 5, 5, 5, 10])

    _Pa, _Qa, _Tna = np.empty(_nx), np.empty(_nx), np.empty(_nx)
    _t_ac, _ = solve_PQ_AC(_a, _b, _c, _d.copy(), np.zeros(_nx), _nx, 21.0, _Pa, _Qa, _Tna)

    _Pf, _Qf, _Tnf = np.empty(_nx), np.empty(_nx), np.empty(_nx)
    _t_fr, _ = solve_PQ(_a, _b, _c, _d.copy(), np.zeros(_nx), _nx, 21.0, 0.0, _Pf, _Qf, _Tnf)

    assert np.allclose(_t_ac, _t_fr), "solve_PQ_AC difiere de solve_PQ"
    tdma_diff = float(np.max(np.abs(_t_ac - _t_fr)))
    return (tdma_diff,)


@app.cell
def _(mo, tdma_diff):
    mo.md(f"""
    ✅ `solve_PQ_AC` ≡ `solve_PQ` — diferencia máxima = `{tdma_diff:g}`
    """)
    return


if __name__ == "__main__":
    app.run()
