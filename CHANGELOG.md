# Changelog

## 0.4.0

- **`joint_web` friendly key in `HollowBlock`**: the half web on the right
  symmetry edge of the cell can now be declared as `joint_web`, mirroring `web`
  on the left edge (the full webs between cavities in the periodic wall are
  `2·web` and `2·joint_web`). Optional — defaults to `web` (uniform webs), so
  results are identical when not declared. The raw key `a12` is still accepted
  and keeps its C-engine meaning, the *full* alternating web:
  `a12 = 2·joint_web`.

## 0.3.0

**Scientific-review release** — every P0 blocker and all P1/P2 items of the
July 2026 documentation/model review are resolved. Results are **not**
bit-reproducible against 0.2.1: the default `hi`, the cavity radiation model,
the 1D layer mapping, the roof-air viscosity, the wall convective constant
and the intermediate-tilt `RF` all changed (each documented below). The
independent validation campaign (EnergyPlus + Borbón hot box) runs in a
separate repository on top of this version.

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

- **Legacy C reference archived**: the original C sources (`legacy_eh/`,
  removed from `main` at "cleaning repo") and the full 0.2.0 development
  history are preserved in the tag `archive/0.2.0-dev`; the C reference will
  live in the validation repository. The four golden tests that need it now
  **skip with a clear message** (and a one-line restore command) instead of
  crashing when `legacy_eh/` is absent; the cached full-day comparison keeps
  running from the committed CSVs.
- **Documentation examples regenerated** with the final physics (radiosity,
  1D layer mapping, air viscosity, Xamán wall constant), default 80×160 mesh:
  hollow-block wall free −0.20 % (`energy_transfer` 47,705 → 47,610), wall AC
  `cooling_energy` −2.26 % (1,611,252 → 1,574,884) and `heating_energy`
  −3.66 %, roof slab free −0.61 % (27,041 → 26,877) J/(m²·day); all converged
  in 2–6 days. `usage.qmd` re-rendered from the new CSVs.
- **Wall-cavity convective constant corrected and de-hardcoded**: the wall
  correlation is the dimensional reduction of Xamán et al. (2005) Eq. (11)
  (turbulent, A = 20, `Nu = 0.0857·Ra^0.3033` — the exponent of `d`,
  `3n−1 = −0.0901`, proves the lineage). The C hardcoded `0.4005`, an
  unrecorded reduction ~0.61× the faithful value; production now computes
  `C_w = 0.0857·k·(gβ/να)^0.3033` at run time from the 300 K property set and
  the configurable air density (≈ 0.589 with the defaults). The legacy 0.4005
  survives only in the C-fidelity golden paths. Wall daily energies shift
  ≈ +1 %; the Borbón hot-box case in the validation campaign tests the
  constant (and the tall-cavity extrapolation) directly.
- **Roof-cavity air viscosity corrected and de-hardcoded**: the C inherited
  `ν = 1.11e-5` m²/s — air at ~240 K, inconsistent with the rest of the
  property set (300 K). ν is now computed as Sutherland `μ(300 K)`
  (1.846e-5 Pa·s, matching Incropera) divided by the **configurable** air
  density, keeping ν and α thermodynamically consistent (`ν ≈ 1.56e-5` with
  the defaults). Only `Slab` + `Fill.AIR` roofs are affected (their Rayleigh
  number was ~43 % overestimated; the Hollands `h_c` decreases slightly).
- **Bibliography completed**: DOIs/issue numbers for Xamán, Hollands, Chow &
  Levermore, Morillón and Borbón (full-text URL); NOM-008-ENER-2001 (Appendix
  B) and NOM-020-ENER-2011 added as standards and cited where the `ho`/`hi`
  defaults are stated; the 3.9 °C (7 °F) long-wave value traced to the
  classical sol–air correction (ASHRAE) with the concept due to Mackey &
  Wright (1944), now cited; the EPW data-dictionary spec cited for the
  hourly-accumulation semantics.
- **New "Assumptions and limits" page**: the canonical assumption →
  consequence table (scope, materials/geometry, surfaces/forcing, cavities),
  linked from the landing page, the Theory menu and both model pages.
- **Average-day factual fixes**: `Ib` is the EPW's **Direct Normal**
  Irradiance (DNI) — it was documented as "beam horizontal"; the GHI/DNI/DHI
  mapping is now stated everywhere, along with the EPW hourly-accumulation
  semantics, the pvlib transposition model and albedo actually used
  (isotropic sky, 0.25 — the defaults), the real `meanDay()` signature and
  the meaning of `day`/`year` on a TMY (synthetic solar date, not a data
  filter).
- **Theory pages homogenized**: `model-1d` and `model-2d` now share the same
  skeleton (Domain and assumptions → Governing equations → Boundary
  conditions → solution modes → shared *average day* → Outputs and units).
  The 2D page gains its explicit assumptions list, a *Solution modes* section
  (width-averaged indoor flux; in AC mode the cavity air keeps floating), an
  outputs table matching the 1D format, and the 1D↔2D axis mapping. The
  incorrect "by periodicity, no heat flows" was replaced by the correct
  argument (the lateral cuts are mirror-symmetry planes, hence adiabatic);
  `Fill.SOLID` is no longer nested under a `Fill.AIR` heading.
- **Linear long-wave factor `RF` (P2)**: `RF` now decreases linearly with the
  surface tilt, `3.9·(1 − tilt/90°)` °C (0 beyond 90°) — the 2016 online
  tool's rule — instead of the binary 3.9-only-at-tilt-0 inherited port.
  Results at `tilt = 0` and `tilt = 90` (all documented cases, and the only
  tilts 2D admits) are bit-identical; intermediate 1D tilts no longer lose
  the sky correction discontinuously.
- **pytz dropped**: `Location.timezone` is now a stdlib `datetime.timezone`
  fixed offset (pytz was only used for `FixedOffset`, and it reached the
  package transitively through pvlib — a dependency pvlib itself is moving
  away from). One less implicit dependency; same offsets, same results.
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
  documented as public API. "Color" is no longer used as a synonym of solar
  absorptance: it remains once, as intuition with a caveat (same-color
  coatings can differ in absorptance). Unit notation unified across docs,
  README and docstrings to the parenthesized SI form: `W/(m·K)`, `W/(m²·K)`,
  `J/(kg·K)`, `J/(m²·day)`.

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
