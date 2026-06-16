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
| — | Extras: `dt` en API · barrido por procesos · `tipo 4` · unidades índices | ⏳ |
| — | Docs: revisar `README.md` (API 2D) · cambio de versión | ⏳ |

> **Lo hecho** vive en las Fases 0–8b (en [`PLAN-2D-hecho.md`](PLAN-2D-hecho.md), cada una con
> su prueba). **Lo pendiente** y sus esquemas/propuestas de integración están agrupados en la
> sección [Trabajo pendiente](#trabajo-pendiente-diseño-y-propuestas) al final.

## Trabajo pendiente (diseño y propuestas)

Lo que falta. El diseño base, los esquemas y la referencia de `System2D` están en
[`PLAN-2D-hecho.md`](PLAN-2D-hecho.md).

### ⏳ Extras opcionales

Refinamientos; ninguno bloquea. (El **motor paralelo `prange`** y la **vigueta en L** a altura
L3+L4+L5+L6 ya están hechos — ver [`PLAN-2D-hecho.md`](PLAN-2D-hecho.md).)

**1. `dt` efectivo en la API** — exponer un paso de tiempo 2D propio en `config2d`/`System2D`
(la Fase 7 midió `dt=10` → ~1.6× con error despreciable) para canjear velocidad↔precisión.

**2. Barrido por procesos** — helper (`solve_many`) o patrón `multiprocessing`/`joblib` para
correr muchas configuraciones independientes en paralelo (cada `solve()` serial); es la palanca
real de velocidad para volumen. ⚠️ No combinar con los hilos de numba.

**3. Variante simétrica / media celda (`tipo 4`)** — optimización geométrica posterior (el resto
del paquete usa celda completa).

**4. Unidades de los índices** — revisión menor pendiente.

**5. Revisar el `README.md`** — integrar la API 2D (`System2D`, `HollowBlock`, `Slab`),
`config2d` (incl. `parallel`), el inspector a escala (`preview`/`section_report`) y el motor
serial/paralelo, con ejemplos de muro y techo.

**6. Cambio de versión** — bump de la versión del paquete que refleje la API 2D (8a/8b) y
actualizar el changelog si aplica.
