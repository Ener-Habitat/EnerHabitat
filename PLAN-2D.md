# PLAN-2D — Pendiente del solver 2D (vigueta y bovedilla)

Portación de `legacy_eh/2dTfree/` al paquete `enerhabitat`. **Este archivo lista solo lo que
falta.** Todo lo **hecho** (Fases 0–8b) y el **diseño/referencia** (mapa del método,
nomenclatura, catálogo de nodos, tolerancias, estructura de archivos, decisiones, riesgos) viven
en [`PLAN-2D-hecho.md`](PLAN-2D-hecho.md).

## Estado del proyecto (hecho ✅ / pendiente ⏳)

| Fase | Qué | Estado |
|------|-----|--------|
| 0 | Golden master del C legacy (standalone + dumps) | ✅ |
| 1 | Geometría/topología (bovedilla rellena) | ✅ |
| 2 | Ensamble de coeficientes `a,b,c,d` | ✅ |
| 3 | Solver de un paso (línea-TDMA + aire interior) | ✅ |
| 4 | Integración temporal + convergencia día-a-día | ✅ |
| 5 | API `System2D` + JIT numba + reducción al 1D | ✅ |
| 6 | Bovedilla con cámara de aire (`tipo 1`) | ✅ |
| 7 | Desempeño y paralelización | ✅ |
| 8a | API producción: **bloque hueco de concreto (muros)** | ✅ |
| 8b | API producción: **vigueta y bovedilla (techos)**, N cavidades, vigueta en L | ✅ |
| 9 | **Aire acondicionado (AC)** en muros y techos 2D (`solveAC`, espejo del 1D) | ✅ |
| — | **Inspector a escala** de la asignación de materiales | ✅ |
| — | Motor serial/**paralelo** (`config2d.parallel`, default serial) · multi-hueco | ✅ |
| — | Rename API español→inglés: `Bovedilla`/`AIRE`/`RELLENA` → `Fill`/`AIR`/`SOLID` | ⏳ |
| — | Extras: `dt` en API · barrido por procesos · `tipo 4` · unidades índices | ⏳ |
| — | `README.md` integra la API 2D (System2D, HollowBlock, Slab, config2d, AC) | ✅ |
| — | Docs: cambio de versión del paquete | ⏳ |

> **Lo hecho** vive en las Fases 0–8b (en [`PLAN-2D-hecho.md`](PLAN-2D-hecho.md), cada una con
> su prueba). **Lo pendiente** y sus esquemas/propuestas de integración están agrupados en la
> sección [Trabajo pendiente](#trabajo-pendiente-diseño-y-propuestas) al final.

## Trabajo pendiente (diseño y propuestas)

Lo que falta. El diseño base, los esquemas y la referencia de `System2D` están en
[`PLAN-2D-hecho.md`](PLAN-2D-hecho.md).

### ⏳ Rename de la API en español → inglés (`Bovedilla`/`AIRE`/`RELLENA`)

**Motivo:** el paquete se publica en inglés, pero el enum **`Bovedilla`**, sus miembros
**`AIRE`/`RELLENA`/`RELLENA_SIMETRICA`**, sus **valores string** (`"aire"/"rellena"/
"rellena_sim"`) y el **parámetro/atributo `bovedilla`** (de `HollowBlock`, `Slab`, `Section2D`,
`SlabSection`) siguen en español. Es la última pieza de API pública en español tras
`colado→topping` y `propiedades/tuplas/__capas`.

**Nombres propuestos (a CONFIRMAR antes de ejecutar):**

| actual | propuesto | nota |
|--------|-----------|------|
| `Bovedilla` (enum) | **`Fill`** | tipo de relleno de la celda del bloque |
| `Bovedilla.AIRE` | **`Fill.AIR`** | |
| `Bovedilla.RELLENA` | **`Fill.SOLID`** | |
| `Bovedilla.RELLENA_SIMETRICA` | **`Fill.SOLID_SYMMETRIC`** | (tipo 4, aún no portado) |
| valores `"aire"/"rellena"/"rellena_sim"` | `"air"/"solid"/"solid_sym"` | solo afectan `.value`/`signature()`/`info()`; no hay datos persistidos |
| parámetro/atributo `bovedilla=` | **`fill_type=`** | **no** usar `fill` a secas: se confunde con `fill_material`. `fill_type` (tipo) + `fill_material` (material) se leen coherentes |

**Decidido con el usuario:** enum **`Fill`** (`AIR`/`SOLID`/`SOLID_SYMMETRIC`) + parámetro
**`fill_type`** — encaja con el `fill_material` ya existente (`fill_type=Fill.SOLID,
fill_material="EPS"`). Alternativas descartadas: `Cavity`/`cavity_type` (choca con la cota
geométrica `cavity`), `FillType`, `cell_type`.

**Decisiones (resueltas):**
1. Nombre del enum y del parámetro: **`Fill` + `fill_type`** (confirmado con el usuario).
2. ¿Alias retro-compatible `Bovedilla`? **No** — rename limpio, sin alias (como `colado→topping`).
3. ¿Renombrar también los valores string? **Sí** (consistencia; afecta solo `.value` interno).

**Touchpoints (inventario):**
- `src/enerhabitat/eh2d.py` (~46): enum + `TIPO_C` + `Section2D.bovedilla`/`build()` +
  `HollowBlock`/`Slab` (param, attr, `signature()`, validación) + `System2D.solve()/solveAC()`
  (`is Bovedilla.AIRE`/`RELLENA`) + `_build_section`.
- `src/enerhabitat/__init__.py`: export `Bovedilla` → `Fill`.
- `src/enerhabitat/ehtools2d.py`: 1 comentario que menciona RELLENA/AIRE.
- **Pruebas** que usan `Bovedilla`: `test_eh2d_slab.py` (15), `test_eh2d_ac.py` (13),
  `test_eh2d_geometry.py` (8), `test_eh2d_hueca.py` (5), `test_eh2d_hollowblock.py` (4),
  `test_eh2d_coeffs.py`, `test_eh2d_fullday.py`, `test_eh2d_step.py`, `test_eh2d_package.py`,
  `test_eh2d_perf.py`.
- `README.md` (~10).

**Pasos:**
1. Definir `Fill` (miembros + valores) en `eh2d.py`; actualizar `TIPO_C`.
2. Renombrar `bovedilla` → `fill_type` en `HollowBlock`, `Slab`, `Section2D`, `SlabSection`
   (param, attr, docstrings, `signature()`), y la validación elemento↔tipo.
3. Actualizar el ruteo en `System2D.solve()/solveAC()` y `_build_section` (`is Fill.AIR`...).
4. Export en `__init__` (`Fill`).
5. Actualizar el comentario en `ehtools2d.py`, las pruebas y el `README.md`.
6. Correr la suite 2D (slab, ac, hollowblock, hueca, inspect, geometry…) — debe pasar sin
   cambios de comportamiento (rename puro).

**Opcional — pase más amplio (identificadores internos con raíz española):** funciones/vars
no públicas `draw_rellena`/`set_krhoc_rellena`/`draw_hueca`/`set_krhoc_hueca`/`_step_hueca(_par)`/
`solve_day_hueca(_ac/_prod/_par)`/`solve_step_hueca`/`Thueco`/`draw_viguetabovedilla*`. No son
API pública (no se exportan); renombrarlas (p. ej. `hueca→cavity`, `rellena→solid`,
`Thueco→T_cav`) es cosmético y mayor. **Fuera del alcance base**; decidir aparte si se quiere
el paquete 100 % en inglés también por dentro.

### ⏳ Extras opcionales

Refinamientos; ninguno bloquea. (El **motor paralelo `prange`**, la **vigueta en L** a altura
L3+L4+L5+L6, el **AC 2D** (`solveAC`) y la integración del **`README.md`** ya están hechos —
ver [`PLAN-2D-hecho.md`](PLAN-2D-hecho.md).)

**1. `dt` efectivo en la API** — exponer un paso de tiempo 2D propio en `config2d`/`System2D`
(la Fase 7 midió `dt=10` → ~1.6× con error despreciable) para canjear velocidad↔precisión.

**2. Barrido por procesos** — helper (`solve_many`) o patrón `multiprocessing`/`joblib` para
correr muchas configuraciones independientes en paralelo (cada `solve()` serial); es la palanca
real de velocidad para volumen. ⚠️ No combinar con los hilos de numba.

**3. Variante simétrica / media celda (`tipo 4`)** — optimización geométrica posterior (el resto
del paquete usa celda completa).

**4. Unidades de los índices** — revisión menor pendiente.

**5. Cambio de versión** — bump de la versión del paquete que refleje la API 2D (8a/8b/9) y
actualizar el changelog si aplica.
