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
| — | Motor **serial** · multi-hueco (`prange` medido a 200×150 → ~break-even, **removido**) | ✅ |
| — | Rename API español→inglés: `Bovedilla`/`AIRE`/`RELLENA` → `Fill`/`AIR`/`SOLID` | ✅ |
| — | Extras: `dt` en API · barrido por procesos · `tipo 4` · unidades índices | ⏳ |
| — | `README.md` integra la API 2D (System2D, HollowBlock, Slab, config2d, AC) | ✅ |
| — | Docs: cambio de versión del paquete (0.2.1, ver `PLAN-README.md`) | ✅ |

> **Lo hecho** vive en las Fases 0–8b (en [`PLAN-2D-hecho.md`](PLAN-2D-hecho.md), cada una con
> su prueba). **Lo pendiente** y sus esquemas/propuestas de integración están agrupados en la
> sección [Trabajo pendiente](#trabajo-pendiente-diseño-y-propuestas) al final.

## Trabajo pendiente (diseño y propuestas)

Lo que falta. El diseño base, los esquemas y la referencia de `System2D` están en
[`PLAN-2D-hecho.md`](PLAN-2D-hecho.md).

### ⏳ Extras opcionales

Refinamientos; ninguno bloquea. (La **vigueta en L** a altura L3+L4+L5+L6, el **AC 2D**
(`solveAC`), el **rename `Bovedilla→Fill`/`bovedilla→fill_type`** y la integración del
**`README.md`** ya están hechos. El **`prange`** se midió y se **removió** por no aportar —
ver [`PLAN-2D-hecho.md`](PLAN-2D-hecho.md).)

**1. `dt` efectivo en la API** — exponer un paso de tiempo 2D propio. **Bajo valor (medido):** el
estudio de `dt` (ver `PLAN-2D-hecho.md`) mostró que `dt` **no afecta el resultado** (1D y 2D) y
que subir `dt`>~10 **ni siquiera acelera** (crecen las iteraciones internas). Casi no aporta.

**2. Barrido por procesos** — helper (`solve_many`) o patrón `multiprocessing`/`joblib` para
correr muchas configuraciones independientes en paralelo (cada `solve()` serial). **Es la palanca
real (medido):** ~6× a 8 procesos, ~10× a 16, escala casi lineal (ver el mapa en
`PLAN-2D-hecho.md`); muy superior al `prange` que se midió (~1.06×) y se removió.

**3. Variante simétrica / media celda (`tipo 4`)** — optimización geométrica posterior (el resto
del paquete usa celda completa).

**4. Unidades de los índices** — revisión menor pendiente.

**5. Cambio de versión** — ✅ hecho: la API 2D salió en **0.2.0** y la documentación
(sitio Quarto, README, CITATION.cff) en **0.2.1** — ver [`PLAN-README.md`](PLAN-README.md).
