# Changelog

## 0.2.0

**2D constructive systems** — cross-sections that are heterogeneous across their width
(concrete hollow-block walls and joist-and-block roofs, with air cavities). The 1D `System`
API is unchanged.

- **`System2D`** — mirror API of `System`: same `Location`, `tilt`, `azimuth`, `absortance`,
  EPW + pvlib `Tsa()`, free-running `solve()` and air-conditioned `solveAC()`. Returns the
  indoor temperature `Ti` and reports `energy_transfer` / `cooling_energy` / `heating_energy`
  (per unit interior area).
- **`HollowBlock`** — concrete hollow block for **walls** (`tilt=90`). The cell can be an
  **air cavity** (`fill_type=Fill.AIR`: wall Nusselt convection + radiation between the cavity
  walls) or **solid-filled** (`fill_type=Fill.SOLID` + `fill_material`, e.g. an insulating core).
- **`Slab`** — joist-and-block roof (*vigueta y bovedilla*) for **roofs** (`tilt=0`): N equal
  cavities, an L-shaped concrete rib, and three solid materials (rib / filler block /
  compression topping); cavities can be air (roof Rayleigh Nusselt) or a solid fill.
- **`solveAC()` for 2D** — cooling/heating demand to hold the indoor air at a setpoint.
- **To-scale section inspector** — `preview()` (matplotlib, optional extra `enerhabitat[viz]`),
  `section_report()`, `section()`.
- **`config2d`** — 2D-only mesh and convergence parameters (`nx`, `ny`, tolerances, `max_days`).
- Serial line-by-line solver (numba JIT), validated against the legacy C golden masters and
  reducing to the 1D path for a homogeneous layer.
- Documentation in English; README gains a 2D section with to-scale figures.

## 0.1.x

1D multilayer opaque constructive systems: `Location`, `System` (free-running and
air-conditioned), global `config`, EPW + pvlib sun-air temperature, finite-volume TDMA solver
over an average day.
