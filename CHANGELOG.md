# Changelog

## Unreleased

- **`Fill.SOLID_SYMMETRIC` removed (P1)**: the C's `tipo 4` (symmetric half
  cell) is deliberately not ported — it is redundant with `Fill.SOLID` on the
  full cell and inapplicable to `Slab` (L-shaped rib breaks the symmetry). The
  member was never functional (it always raised `NotImplementedError`), so no
  working code is affected; the docs always listed two values.
- **Targeted range validation (P1)**: physically corrupting inputs now raise
  `ValueError` instead of being accepted silently — `absortance` outside
  `[0, 1]` (1D and 2D), `config.La/ho/hi`/air properties ≤ 0, non-integer or
  < 3 `Nx`, layer thicknesses ≤ 0, materials with `k/rho/c` ≤ 0 (rejected at
  `.ini` load, keeping the previous materials), and `config2d.validate()`
  (mesh/tolerances/caps) run at the start of every 2D solve. Missing layer
  materials now raise a `KeyError` that lists the available names.
- **Explicit errors for the materials file (P1)**: assigning a missing path to
  `config.file` now raises `FileNotFoundError` (it used to print and continue,
  leaving internal attributes unset and crashing later with a cryptic
  `AttributeError`); the previously loaded materials are kept on failure.
  Importing EnerHabitat without a `materials.ini` in the working directory is
  now silent and leaves `config.materials == {}`; the `file` getter no longer
  performs I/O nor returns `None`.

- **Fractional UTC offsets fixed (P2)**: the EPW header's decimal UTC offset
  is now preserved via `pytz.FixedOffset` — it used to be truncated
  (`+5.5 → +5`, up to 45 min of solar-time error in half/quarter-hour zones
  such as India, Nepal or Newfoundland) and encoded as `Etc/GMT±N`, which
  cannot represent fractions at all. Integer-offset EPWs (e.g. Mexico) are
  bit-identical. `Location.timezone` is documented as a fixed-offset zone
  (local standard time, no DST).
- **API-accuracy doc fixes (P1 series)**: `solve()`/`solveAC()` return a
  `pandas.Series` named `"Ti"`, not a DataFrame — corrected in usage/api/
  model-1d pages, 1D docstrings and type annotations (which also promised a
  nonexistent `energy` argument and `Qcool, Qheat` return values; the energies
  are instance attributes). Use `.to_frame("Ti")` when a DataFrame is needed.
  The `Tsa()` cache description was unified across pages (setters invalidate,
  the next `Tsa()`/solver call recomputes; a prior manual call is optional and
  `layers` does not affect `Tsa`), with a new warning: mutating `layers` in
  place bypasses the **solver** cache invalidation — assign the list or use
  `add_layer()`/`remove_layer()`. `config.reset()` is now documented as
  restoring the *numeric* defaults only: `file`/`materials` are user-provided
  and are kept. The `System2D` ↔ `System` differences are now spelled out
  (2D-only `setpoint`, missing `flag()`, `copy()` shares the `Location` in 2D
  but re-creates it in 1D, writable vs read-only result attributes, per-mode
  caches, extra 2D outputs); README no longer claims the interfaces are
  identical. The 2D 7-layer limit (fixed `L1…L7` slots inherited from the C
  engine; already enforced with a clear `ValueError`) is now documented in
  api/usage. The 2D outputs `solve_dataframe` (columns, units, `Thueco`
  semantics per element, and the C-inherited sampling convention: `Tso` at the
  start of each step, `Ti`/`Tsi` at the end), `Tfield` and `Qout` are now
  documented as public API.

- **Cavity radiation solved with radiosity** (P0-02): the grey diffuse
  enclosure of each cavity is now solved exactly via Gebhart transfer factors
  `𝔉 = ε²(I−(1−ε)F)⁻¹F`, precomputed once per geometry and fed to the
  unchanged kernels (which evaluate the same pairwise form with `E = 1`). The
  previous direct-exchange approximation `ε·σ·F·ΔT⁴` — inherited from the C —
  overestimated the dominant pair by ~10 % at ε = 0.9. **2D results with
  `Fill.AIR` change**; `ε = 1` and the legacy/golden paths are bit-identical.
  Per-surface emissivities are now supported by the formulation.

- **1D layer-to-mesh mapping corrected** (P0-05): layers are now assigned by
  cumulative coordinates (cell-centre material, thickness-weighted ρc — total
  thermal mass conserved exactly) and interior face conductances are computed
  from the exact series resistance between cell centres, so material
  interfaces may fall anywhere inside a cell. The previous `int(L/dx)+1`
  counting gave each layer one extra cell (e.g. two 50 mm layers at Nx=10
  became 60/40 mm), made results depend on layer ordering, and either deleted
  a sub-cell sheet placed first (a 0.5 mm metal sheet vanished silently) or
  turned it into a fake material. **Multilayer 1D results change** (toward the
  correct geometry); single-layer systems are bit-identical.

- **Daily (periodic) convergence over all states** (P0-04): the day-to-day
  criterion now closes **every persisted state** — the solid field, the indoor
  air `T_i` (free-running) and each cavity air `T_h` — via
  `day_error = max(C_solid, |ΔT_i|, max|ΔT_h|) ≤ tol_day`. The 1D solver gets a
  `MAX_DAYS = 60` cap (its loop previously had none) and `System` (1D) now
  exposes `days`, `day_error`, `converged` and `energy_imbalance`, warning with
  a `RuntimeWarning` when not converged, like `System2D`.
- **Energy-closure diagnostic**: free-running solves report
  `energy_imbalance = |Qin − Qout| / max(Qin, Qout)` (≈ 0 in the periodic
  regime) in both 1D and 2D.
- **Docs**: the claim that the AC mode is "purely implicit" was corrected —
  with `Fill.AIR` the cavity air still advances explicitly in AC mode; its
  stability parameter `λ_h·Δt = h_c·P_cav·Δt/(ρ_a·c_a·A_cav)` is now documented
  (≈ 0.2 wall cavity, up to ≈ 1 roof cavities at Δt = 10 s).

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
