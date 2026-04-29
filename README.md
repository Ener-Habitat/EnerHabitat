# EnerHabitat

[![PyPI version](https://img.shields.io/pypi/v/enerhabitat.svg)](https://pypi.org/project/enerhabitat/)
[![Python versions](https://img.shields.io/pypi/pyversions/enerhabitat.svg)](https://pypi.org/project/enerhabitat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/Ener-Habitat/EnerHabitat/blob/main/LICENSE)

**EnerHabitat** is a Python package for the thermal simulation of opaque
constructive systems (walls and roofs) driven by EPW weather data. It solves
the one-dimensional, time-dependent heat conduction equation across multi-layer
systems and produces hourly indoor temperatures and air-conditioning energy
demands for an *average day* of a chosen month.

## Contents
- [Overview](#overview)
- [Theoretical background](#theoretical-background)
- [Installation](#installation)
- [Recommended folder structure](#recommended-folder-structure)
- [Key concepts](#key-concepts)
- [Quickstart](#quickstart)
- [Workflow](#workflow)
- [Examples](#examples)
  - [Two-layer system without air conditioning](#two-layer-system-without-air-conditioning)
  - [Two-layer system with air conditioning](#two-layer-system-with-air-conditioning)
- [API reference](#api-reference)
  - [Location](#location)
  - [System](#system)
- [Config (global)](#config-global)
- [Materials](#materials)
- [Dependencies](#dependencies)
- [Authors](#authors)
- [License](#license)

## Overview

EnerHabitat models the heat transfer through opaque constructive systems
**without windows, ventilation or infiltration**. Each layer of the system is
described by a material name and three thermal properties:

- thermal conductivity `k` (W/m·K)
- density `rho` (kg/m³)
- specific heat `c` (J/kg·K) — written `c_p` in the heat equation below

These three names (`k`, `rho`, `c`) are the **exact keys** expected in
`materials.ini`; they are case-sensitive and Greek letters are not accepted.

Given an EPW file and a constructive system, EnerHabitat computes:

| Symbol | Description |
| ------ | ----------- |
| `Ta`   | Outdoor ambient temperature |
| `Tsa`  | Sun-air temperature |
| `Ti`   | Indoor temperature |
| `Tn`   | Adaptive comfort (neutral) temperature |
| `Ig`   | Global horizontal irradiance |
| `Ib`   | Direct normal irradiance |
| `Id`   | Diffuse horizontal irradiance |
| `Is`   | Solar irradiance on the tilted surface |

## Theoretical background

EnerHabitat solves the 1-D, time-dependent heat conduction equation across the
constructive system:

```
              ∂T          ∂²T
    ρ   c_p   ──  =  k   ────
              ∂t          ∂x²
```

The exterior boundary condition uses the **sun-air temperature**, which combines
convection, short-wave solar gain and long-wave radiative losses:

```
    T_sa = T_o + (I_s   a) / h_o + RF
```

where:

- `T_o` — outdoor ambient temperature
- `I_s` — solar irradiance incident on the surface
- `a`   — external solar absorptance
- `h_o` — outdoor convective heat transfer coefficient
- `RF`  — long-wave radiative loss factor

The equation is discretised with **finite control volumes** and solved with the
**TDMA** (Tri-Diagonal Matrix Algorithm). The simulation runs over an
**average day** of a selected month — built from the EPW — and is iterated
until a periodic (oscillatory) steady state is reached.

Two solution modes are available:

- **Free-running** — `solve()`: no air conditioning is applied; the indoor
  temperature follows the dynamics of the constructive system.
- **Air-conditioned** — `solveAC()`: the indoor temperature is held at a
  comfort setpoint derived from the **Humphreys & Nicol** adaptive comfort
  model combined with **Morillón's** comfort-zone amplitude proposal.
  EnerHabitat then applies the cooling or heating needed at every time step to
  keep `Ti` at that setpoint, and reports the resulting `cooling_energy` and
  `heating_energy` demands.

## Installation

```bash
pip install enerhabitat
```

With [uv](https://docs.astral.sh/uv/) (we love it and warmly encourage its use
— it is fast, reproducible, and our recommended way to install EnerHabitat):

```bash
uv add enerhabitat
```

EnerHabitat requires **Python ≥ 3.10**.

## Recommended folder structure

`materials.ini` is **required** — EnerHabitat ships with no default materials,
so you must provide this file (see [Materials](#materials) for its format).

```
project/
├── main.py
├── materials.ini      # Material properties (REQUIRED — user-provided)
└── epw/
    ├── ...
    └── example.epw
```

## Key concepts

- **`Location`** reads an EPW file and computes the average day with `meanDay()`.
- **`System`** combines a `Location` and a list of layers and computes `Tsa()`,
  `solve()` and `solveAC()`.
- **`config`** is a global instance whose attributes (materials file,
  discretisation, convection coefficients, time step) affect every subsequent
  computation.

## Quickstart

EnerHabitat does **not** ship with pre-loaded materials. Before running anything,
create a `materials.ini` file (see [Materials](#materials)) in your working
directory or point `eh.config.file` to its location.

```python
import enerhabitat as eh

# 1) Materials file (required — no defaults are bundled)
eh.config.file = "./materials.ini"

# 2) Location from an EPW file
loc = eh.Location("./epw/example.epw")

# 3) Define the constructive system
wall = eh.System(location=loc)
wall.azimuth = 90
wall.absortance = 0.3
wall.layers = [("Adobe", 0.20)]      # outside → inside

# 4) Average day and solar inputs
loc.meanDay(month=5, year=2025)
wall.Tsa()

# 5) Solve
ti = wall.solve()
```

## Workflow

To simulate a wall (or roof) you need to:

1. **Geolocate it** — pass an EPW file to `Location`.
2. **Orient it** — set `azimuth` (and `tilt` if needed).
3. **Define its color** — set `absortance`.
4. **Define its layers** — set `layers` from outside to inside.
5. **Choose the period** — call `location.meanDay(month, year)`.
6. **Compute `Tsa()`**, then choose one solver:
   - **`solve()`** — *without* air conditioning (free-running): the indoor
     temperature `Ti` evolves freely with the dynamics of the constructive
     system.
   - **`solveAC()`** — *with* air conditioning: the indoor temperature is held
     at a comfort setpoint and the cooling/heating energy required is
     reported.

Both `solve()` and `solveAC()` return pandas DataFrames indexed by time of day.

## Examples

### Two-layer system without air conditioning

```python
import enerhabitat as eh
import pandas as pd

epw_file = "epw/MEX_CAM_Campeche-Ignacio.766961_TMYx.epw"

wall = eh.System(eh.Location(epw_file))
wall.azimuth = 90
wall.absortance = 0.3
wall.layers = [("Mortero", 0.025), ("Ladrillo", 0.10)]
wall.location.meanDay(month=5, year=2025)
wall.Tsa()

# Free-running solution
data = wall.solve()

# Optionally attach Tsa to the result (useful when color and location
# do not change between runs)
data = pd.concat([data, wall.Tsa().asfreq("10min")], axis=1)
```

### Two-layer system with air conditioning

```python
import enerhabitat as eh
import pandas as pd

epw_file = "epw/MEX_CAM_Campeche-Ignacio.766961_TMYx.epw"

wall = eh.System(eh.Location(epw_file))
wall.azimuth = 90
wall.absortance = 0.3
wall.layers = [("Mortero", 0.025), ("Ladrillo", 0.10)]
wall.location.meanDay(month=5, year=2025)
wall.Tsa()

# Air-conditioned solution: setpoint at the upper comfort bound
data = wall.solveAC()
data = pd.concat([data, wall.Tsa().asfreq("10min")], axis=1)

# Cooling and heating energy demands, in J/(m²·day) over one average day
print(wall.cooling_energy, wall.heating_energy)
```

## API reference

```python
import enerhabitat as eh
```

### Location

```python
loc = eh.Location("./epw/example.epw")
```

**Attributes**

The EPW path is stored in `file`. The following are read-only and recovered
from the EPW header — change `file` to update them:

- `city` — `str`, city from the EPW header
- `latitude` — `float`, degrees
- `longitude` — `float`, degrees
- `altitude` — `float`, metres
- `timezone` — `pytz.timezone`

```python
loc.file = "./epw/other.epw"
```

**Methods**

- `meanDay(month, year)` — average-day DataFrame (`Ta`, `Ig`, `Ib`, `Id`, `Tn`)
- `copy()` — returns a copy of the instance
- `info()` — prints instance attributes
- `flag()` — `dict` with metadata of the last `meanDay()` call

```python
loc.meanDay(month=6, year=2020).info()
loc.info()
print(loc.flag()["date"])
```

### System

```python
loc = eh.Location("./epw/example.epw")
wall = eh.System(location=loc)
```

**Attributes**

- `location` — associated `Location`
- `tilt` — `float`, degrees from horizontal (`0` = roof, `90` = vertical wall)
- `azimuth` — `float`, surface azimuth in degrees (pvlib convention, clockwise
  from north):

  | Direction | Azimuth |
  | --------- | ------- |
  | North     | `0`     |
  | East      | `90`    |
  | South     | `180`   |
  | West      | `270`   |
- `absortance` — `float` in `[0, 1]`
- `layers` — `list[tuple[str, float]]` of `(material, thickness_m)`,
  ordered **from outside to inside**

```python
wall.location = loc_2
wall.tilt = 0
wall.azimuth = 45
wall.absortance = 0.3

wall.layers = [("Adobe", 0.10), ("Acero", 0.05), ("Ladrillo", 0.02)]

wall.add_layer("Mortero", 0.20)   # appended at the inside
wall.remove_layer(2)              # removes layer at index 2
```

Read-only result attributes (all expressed in **J/(m²·day)** — energy per
unit surface area, accumulated over one converged average day):

- `energy_transfer` — total energy transferred to the indoor side from
  `solve()`
- `heating_energy` — heating demand from `solveAC()`
- `cooling_energy` — cooling demand from `solveAC()`

> Units: `hi · Δt · ΔT` with `hi` in W/(m²·K), `Δt` in seconds and `ΔT` in K
> yields **J/m²**, and the loop accumulates these contributions over the 24 h
> of the average day, so the reported value is **J/(m²·day)**.
> Divide by `3600` to get Wh/(m²·day) or by `3.6e6` to get kWh/(m²·day).

**Methods**

- `Tsa()` — sun-air temperature and `Is` from `Location.meanDay()`
- `solve()` — indoor temperature `Ti` (free-running)
- `solveAC()` — cooling and heating energy with constant indoor setpoint
- `copy()` — returns a copy of the instance
- `info()` — prints attributes
- `flag()` — reports whether the cached value was recomputed

```python
wall.Tsa().info()

ti = wall.solve()
energy = wall.energy_transfer

wall.solveAC()
c_energy = wall.cooling_energy
h_energy = wall.heating_energy
```

> Tip: when color and location do not change between runs, attach the sun-air
> temperature to the result with
> `data = pd.concat([data, wall.Tsa().asfreq("10min")], axis=1)`.

## Config (global)

`config` is a global singleton that stores parameters shared by every
`Location` and `System`. Changing it affects **all** subsequent computations.

**Attributes**

- `file` — path to the `.ini` file with material properties
- `La` — length of the fictional indoor space (m)
- `Nx` — number of control volumes used to discretise the system
- `ho` — outdoor convective heat transfer coefficient (W/m²·K)
- `hi` — indoor convective heat transfer coefficient (W/m²·K)
- `dt` — time step (seconds)

```python
eh.config.file = "./materials.ini"

eh.config.La = 2.0
eh.config.Nx = 300
eh.config.ho = 12
eh.config.hi = 8.3
eh.config.dt = 60
```

`config.materials` is a read-only `dict` keyed by material name:

```python
adobe = eh.config.materials["Adobe"]
adobe.k     # W/m·K
adobe.rho   # kg/m³
adobe.c     # J/kg·K
```

**Methods**

- `info()` — prints current `config` values
- `to_dict()` — returns parameters as a `dict`
- `reset()` — restores default values
- `materials_list()` — list of material names defined in `file`
- `materials_dict()` — dict of material properties

## Materials

EnerHabitat **does not bundle any default materials**. You must supply a
`materials.ini` file — by default the package looks for `materials.ini` in the
current working directory; otherwise set `eh.config.file` to the path you want
to use.

Material properties are declared in an `.ini` file, with the material name as
the section header and `k`, `rho` and `c` as keys:

```ini
[concrete]
k   = 1.35   # Thermal conductivity, W/m·K
rho = 1800   # Density, kg/m³
c   = 1000   # Specific heat, J/kg·K

[adobe]
k   = 0.58
rho = 1500
c   = 1480
```

Point `config.file` to a different file when you need to switch material sets:

```python
eh.config.file = "./config/new_materials.ini"
```

If `config.file` points at a missing file, EnerHabitat will report
`Error: <path> not found` and `materials` will be empty — `System.solve()` will
fail because the layer materials cannot be resolved.

## Dependencies

Direct dependencies (declared in `pyproject.toml`):

- [numba](https://numba.pydata.org/)
- [pvlib](https://pvlib-python.readthedocs.io/)

Pulled in transitively and used internally:

- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [pytz](https://pypi.org/project/pytz/)

## Authors

Developed at the **Instituto de Energías Renovables, UNAM**.

- Guillermo Barrios del Valle — <gbv@ier.unam.mx>
- Fernando Rodríguez Calderón — <ferrodriguez2509@gmail.com>

Source code: <https://github.com/Ener-Habitat/EnerHabitat>
Issues: <https://github.com/Ener-Habitat/EnerHabitat/issues>

## License

Released under the [MIT License](https://github.com/Ener-Habitat/EnerHabitat/blob/main/LICENSE).
