# EnerHabitat

[![PyPI version](https://img.shields.io/pypi/v/enerhabitat.svg)](https://pypi.org/project/enerhabitat/)
[![Python versions](https://img.shields.io/pypi/pyversions/enerhabitat.svg)](https://pypi.org/project/enerhabitat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/Ener-Habitat/EnerHabitat/blob/main/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-ener--habitat.github.io-blue)](https://ener-habitat.github.io/EnerHabitat/)

**EnerHabitat** is a Python package for the thermal simulation of opaque
constructive systems (walls and roofs) driven by EPW weather data. It solves
the time-dependent heat conduction equation across multi-layer systems — in
**1D** for homogeneous layers (`System`) and in **2D** for units that are
heterogeneous across their width, such as concrete hollow-block walls and
joist-and-block (*vigueta y bovedilla*) roofs (`System2D`) — and produces
indoor temperatures and air-conditioning energy demands for an *average day*
of a chosen month.

## 📖 Documentation

Full documentation lives at **<https://ener-habitat.github.io/EnerHabitat/>**:

- [Usage](https://ener-habitat.github.io/EnerHabitat/usage.html) — workflow,
  complete examples (1D/2D, free-running and air-conditioned), configuration
  and materials.
- Theory —
  [1D model](https://ener-habitat.github.io/EnerHabitat/model-1d.html) ·
  [2D model](https://ener-habitat.github.io/EnerHabitat/model-2d.html) ·
  [Numerical method](https://ener-habitat.github.io/EnerHabitat/numerics.html)
  (equations, boundary conditions, convergence, validation).
- [API reference](https://ener-habitat.github.io/EnerHabitat/api.html).

## Overview

EnerHabitat models the heat transfer through opaque constructive systems
**without windows, ventilation, infiltration or internal heat gains**. Each
layer is described by a material name and three properties: thermal
conductivity `k` (W/(m·K)), density `rho` (kg/m³) and specific heat `c`
(J/(kg·K)) — supplied by a user `materials.ini` file (no defaults are bundled).

Given an EPW file and a constructive system, EnerHabitat computes the outdoor
(`Ta`), sun–air (`Tsa`), indoor (`Ti`) and neutrality (`Tn`) temperatures,
the solar irradiances (`Ig`, `Ib`, `Id`, `Is`) and the energy demands over
the average day of the selected month.

## Theoretical background (summary)

The temperature field in each layer obeys the 1D time-dependent heat
conduction equation,

$$
\rho\, c_p\, \frac{\partial T}{\partial t} = k\, \frac{\partial^2 T}{\partial x^2},
$$

with flux continuity at layer joints. At the outdoor surface the boundary
condition uses the **sun–air temperature**, which lumps convection, absorbed
solar radiation and the long-wave sky exchange:

$$
T_{sa} = T_a + \frac{a\, I_s}{h_o} - RF,
$$

with `a` the solar absorptance, $I_s$ the irradiance on the tilted surface
(computed with pvlib) and $RF$ decreasing linearly from 3.9 °C at `tilt = 0`
(roof) to 0 at `tilt = 90` (wall).
At the indoor surface the system exchanges heat with the indoor air,
and two solution modes exist:

- **Free-running** — `solve()`: the indoor air is a lumped thermal mass whose
  temperature `Ti` evolves freely; the daily energy delivered to it is
  reported as `energy_transfer`.
- **Air-conditioned** — `solveAC()`: `Ti` is held at the neutrality
  temperature of the adaptive comfort model of Humphreys & Nicol,
  $T_n = 0.54\,\overline{T_a} + 13.5$ °C, and the required `cooling_energy`
  and `heating_energy` are reported. (The average-day data also includes the
  comfort-zone half-width `DeltaTn`, after Morillón, for comfort analyses.)

For 2D systems, `System2D` solves the same problem on the unit's
cross-section, adding the cavity physics: radiation between the cavity walls
and temperature-dependent Nusselt convection with a lumped cavity-air node.

The equations are discretised with implicit **finite control volumes** and
solved with the **TDMA**; the average day is iterated until the solution is
periodic. Full derivations, boundary conditions, convergence criteria and the
validation record are in the
[theory pages](https://ener-habitat.github.io/EnerHabitat/model-1d.html).

## Installation

```bash
pip install enerhabitat
```

With [uv](https://docs.astral.sh/uv/) (we love it and warmly encourage its
use — fast, reproducible, and our recommended way to install EnerHabitat):

```bash
uv add enerhabitat
```

EnerHabitat requires **Python ≥ 3.10**. The section inspector plots are an
optional extra: `pip install enerhabitat[viz]`.

## Quickstart

EnerHabitat ships **no** materials: create a `materials.ini` in your working
directory (or point `eh.config.file` to one) before running anything. A
minimal file for this example:

```ini
[Adobe]
k   = 0.58    # W/(m·K)
rho = 1500    # kg/m³
c   = 1480    # J/(kg·K)
```

```python
import enerhabitat as eh

# 1) Materials file (required — no defaults are bundled)
eh.config.file = "./materials.ini"

# 2) Location from an EPW file
loc = eh.Location("./epw/example.epw")

# 3) Define the constructive system
wall = eh.System(location=loc)
wall.azimuth = 90                    # east-facing
wall.absortance = 0.3
wall.layers = [("Adobe", 0.20)]      # outside → inside

# 4) Average day and solar inputs
loc.meanDay(month=5, year=2025)
wall.Tsa()

# 5) Solve (free-running); Tsa() and solve() share the same time grid,
#    so results concatenate directly.
ti = wall.solve()
print(wall.energy_transfer)          # J/(m²·day)
```

For a wall with air conditioning, call `wall.solveAC()` and read
`wall.cooling_energy` / `wall.heating_energy`.

### 2D systems

`System2D` is used like `System` (see the
[API page](https://ener-habitat.github.io/EnerHabitat/api.html#system2d) for the
differences); its `layers` list contains
**exactly one** 2D element — a `HollowBlock` (walls, `tilt = 90`) or a `Slab`
(joist-and-block roofs, `tilt = 0`). The materials named below (`Concreto`,
`Mortero`, `Yeso`) must also be defined in your `materials.ini` (see the
[full example set](https://ener-habitat.github.io/EnerHabitat/usage.html#materials-file)):

```python
block = eh.HollowBlock(
    material   = "Concreto",
    emissivity = 0.9,
    geometry   = {"web": 0.02, "block_width": 0.16,
                  "cover_top": 0.02, "cavity": 0.08, "cover_bottom": 0.02},
)

wall = eh.System2D(eh.Location("./epw/example.epw"))
wall.tilt = 90
wall.azimuth = 90
wall.absortance = 0.6
wall.layers = [("Mortero", 0.02), block, ("Yeso", 0.01)]
wall.location.meanDay(month=5, year=2025)
wall.Tsa()
ti = wall.solve()
```

All the examples — the full 1D/2D × free-running/AC matrix, the
joist-and-block roof, and the to-scale section inspector — are in the
[Usage page](https://ener-habitat.github.io/EnerHabitat/usage.html).

## API at a glance

| Object | Purpose | Key methods / attributes |
| ------ | ------- | ------------------------ |
| `Location` | Reads an EPW file, builds the average day | `meanDay(month, year)` |
| `System` | 1D multilayer wall/roof | `layers`, `Tsa()`, `solve()`, `solveAC()`, `energy_transfer`, `cooling_energy`, `heating_energy` |
| `System2D` | 2D heterogeneous wall/roof | same interface as `System`, plus `preview()`, `section_report()`, `days` |
| `HollowBlock` / `Slab` | The 2D element inside `System2D.layers` | `material(s)`, `fill_type` (`Fill.AIR`/`Fill.SOLID`), `geometry` |
| `config` | Global parameters | `file`, `La`, `Nx`, `ho`, `hi`, `dt` *(fixed)* |
| `config2d` | 2D mesh & convergence | `nx`, `ny`, `tol_inner`, `tol_day`, `max_days` |

Defaults for `ho` (13) and `hi` (8.1 W/(m²·K)) are the NOM-008/020-ENER values
(`hi` is the vertical-surface value, applied to all orientations);
`dt` is fixed at 10 s (see
[why](https://ener-habitat.github.io/EnerHabitat/numerics.html#indoor-air-coupling-free-running-mode)).
Full reference:
[API page](https://ener-habitat.github.io/EnerHabitat/api.html).

## Dependencies

Direct dependencies: [numba](https://numba.pydata.org/) and
[pvlib](https://pvlib-python.readthedocs.io/) (numpy, pandas and pytz come
with them). Optional: matplotlib via `enerhabitat[viz]`.

## How to cite

If you use EnerHabitat in academic work, please cite the reference paper:

> Barrios, G., Casas, J.M., Huelsz, G., Rojas, J. (2016). *Ener-Habitat: An
> online numerical tool to evaluate the thermal performance of homogeneous
> and non-homogeneous envelope walls/roofs*. Solar Energy 131, 296–304.
> <https://doi.org/10.1016/j.solener.2015.12.017>

```bibtex
@article{Barrios2016,
  author  = {Barrios, G. and Casas, J.M. and Huelsz, G. and Rojas, J.},
  title   = {Ener-Habitat: An online numerical tool to evaluate the thermal
             performance of homogeneous and non-homogeneous envelope walls/roofs},
  journal = {Solar Energy},
  volume  = {131},
  pages   = {296--304},
  year    = {2016},
  doi     = {10.1016/j.solener.2015.12.017}
}
```

The repository also ships a [`CITATION.cff`](https://github.com/Ener-Habitat/EnerHabitat/blob/main/CITATION.cff)
(GitHub's *Cite this repository* button) covering the software itself.

## Authors

Developed at the **Instituto de Energías Renovables, UNAM**.

- Guillermo Barrios del Valle — <gbv@ier.unam.mx>
- Fernando Rodríguez Calderón — <ferrodriguez2509@gmail.com>

Source code: <https://github.com/Ener-Habitat/EnerHabitat> ·
Issues: <https://github.com/Ener-Habitat/EnerHabitat/issues>

## License

Released under the [MIT License](https://github.com/Ener-Habitat/EnerHabitat/blob/main/LICENSE).
