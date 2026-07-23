# Changelog

## Unreleased

- **2D inner convergence criterion** (P0-03 of the documentation review): the
  production solvers replace the C-inherited *signed* mean relative change —
  which can stop early by cancellation and divides by temperatures in °C —
  with two non-cancellable checks in °C: the max node update of the last sweep
  **and** the max scaled residual of the discrete equations
  (`|a_P·T_P − Σa_nb·T_nb − b| / a_P`), the residual-based monitor recommended
  by Patankar (1980). `config2d.tol_inner` now means this tolerance
  (default `1e-8` °C).
- **Convergence diagnostics**: new `config2d.max_inner` cap (10⁴ sweeps/step);
  `System2D` exposes `converged`, `day_error` and `inner_iterations` after each
  solve and emits a `RuntimeWarning` when a solve does not converge.
- The C-faithful ports and golden-master tests keep the legacy criterion
  (`legacy=True`) and still reproduce the C bit-for-bit (field and iteration
  counts).
- Default indoor film coefficient corrected to the NOM vertical-surface value:
  `config.hi = 8.1` W/(m²·K) (was 8.6, wrongly attributed to the NOM).

## 0.2.1

**Documentation release** — no changes to the solvers or the API.

- **Documentation site**: <https://ener-habitat.github.io/EnerHabitat/> (Quarto on
  GitHub Pages). Theory pages with numbered equations and bibliography — 1D model
  (physical problem, boundary conditions), 2D model (cavity physics: radiation, Nusselt
  correlations, lumped cavity air) and numerical method (implicit control volumes,
  TDMA, convergence criteria, differences with the 2016 *Solar Energy* paper,
  validation record).
- **Executable examples**: the 1D examples run at render time; the 2D ones (hollow
  block free/AC, joist-and-block roof) are pre-computed by `docs/run_examples.py` with
  the default mesh and shown with their plots, energies and compute times.
- **README** restructured as a landing page (usage summary + links to the site);
  fixed the sign convention of `RF` in the sun–air temperature and the description of
  the AC setpoint (`Tn`; `DeltaTn` is data only).
- **`CITATION.cff`** (GitHub *Cite this repository*) and a *How to cite* section.
- `Documentation` URL in the PyPI sidebar (`[project.urls]`).

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
