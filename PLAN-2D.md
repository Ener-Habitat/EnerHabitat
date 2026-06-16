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
| 9 | **Aire acondicionado (AC)** en muros y techos 2D (`solveAC`, espejo del 1D) | ⏳ |
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

### ⏳ Fase 9 — Aire acondicionado (AC) en muros y techos 2D

**Objetivo:** dar a `System2D` un `solveAC()` **espejo del `System.solveAC()` 1D**, para el
bloque hueco (muro) y la vigueta y bovedilla (techo) — calcular la **energía de enfriamiento y
calentamiento** necesaria para mantener el aire interior en un setpoint. Reúsa toda la
infraestructura 8a/8b; el cambio por kernel es chico.

**Flujo de trabajo del usuario (idéntico 1D ↔ 2D):** el único cambio es `solve()` → `solveAC()`
y leer `cooling_energy`/`heating_energy` en vez de `energy_transfer`.

```python
import enerhabitat as eh
eh.config.file = "./materials.ini"
loc = eh.Location("./epw/example.epw")
loc.meanDay(month=5, year=2025)

# --- 1D (vigente) ---
sys = eh.System(loc, tilt=90, azimuth=90, absortance=0.6)
sys.layers = [("Concreto", 0.15), ("EPS", 0.05)]
sys.Tsa()
ti = sys.solveAC()                       # mantiene Tint en el setpoint
print(sys.cooling_energy, sys.heating_energy)   # Qcool, Qheat (energy_transfer = None)

# --- 2D (objetivo, Fase 9) — mismo flujo ---
block = eh.HollowBlock("Concreto", emissivity=0.9, geometry={...})
wall = eh.System2D(loc, tilt=90, azimuth=90, absortance=0.6)
wall.layers = [("Aplanado", 0.02), block, ("Yeso", 0.01)]
wall.Tsa()
ti = wall.solveAC()                      # idéntico; Tint fijo en el setpoint
print(wall.cooling_energy, wall.heating_energy)
```

**Flujo interno (qué cambia frente a la flotación libre):**

| paso | flotación libre (`solve`) | aire acondicionado (`solveAC`) |
|------|---------------------------|--------------------------------|
| aire interior `Tint` | nodo lumped, se **integra** cada paso | **fijo** en el setpoint (no se integra) |
| frontera interior `j=ny-1` | `hi·dx·Tint` con `Tint` que evoluciona | `hi·dx·Tset` con `Tset` constante |
| aire del hueco `Thueco` (AIRE) | flota | **flota igual** (el AC no toca la cavidad) |
| energía | `Qin`/`Qout` (en régimen `Qin≈Qout`) | `Qcool` (flujo neto entra) / `Qheat` (sale) |
| salida | `energy_transfer=Qin` | `cooling_energy=Qcool`, `heating_energy=Qheat` |

Lazo interno (ambos): convergencia día-a-día → por paso, ensamblar coeficientes con `Tsa[t]` y
la frontera interior, resolver el campo (línea-TDMA), y — en libre integrar `Tint`; en AC dejar
`Tset` fijo y **acumular** `Qcool/Qheat` según el signo del flujo neto en la superficie interior.

**Modelo (idéntico al 1D, `ehframe.__calc_solve(AC=True)`):**
- El aire interior `Tint` se **fija en un setpoint** (default `Tn.mean()`, como el 1D; opción de
  exponer un atributo `setpoint`/`Tset`). **No** se integra en el tiempo (a diferencia de la
  flotación libre, donde `Tint` es un nodo lumped que se marcha paso a paso).
- Se resuelve el campo 2D en cada paso con la frontera interior usando ese `Tint` fijo. En
  **AIRE** (hueca/slab), el aire del hueco `Thueco` (o `Thueco[c]` por cavidad) **sigue
  flotando** — el AC solo controla el aire del **recinto** (`Tint`), no la cavidad.
- **Carga del AC por paso** = flujo neto en la superficie interior, con el mismo criterio de
  signo del 1D: `e = (Σ_i hi·dx·(T[i,ny-1] − Tint))·dt / X`; si `e>0` → enfriamiento
  (`Qcool += e`), si `e<0` → calentamiento (`Qheat -= e`). (El 1D parte el signo sobre su único
  nodo de superficie; en 2D es el **flujo neto** sobre los `nx` nodos interiores, que reduce al
  1D para capa homogénea.)
- **Convergencia día-a-día** (`C > tol_day`) hasta permanente, igual que la flotación libre.

**Kernels (`ehtools2d`), variante `_ac` de cada motor** (serial + paralelo, respetando
`config2d.parallel`):
- `solve_day_2d_ac` — RELLENA / conducción pura.
- `solve_day_hueca_ac` — muro (`Thueco` flota, `Tint` fijo).
- `solve_day_slab_ac` — techo N cavidades (`Thueco[c]` flotan, `Tint` fijo).
Cada uno = su versión de producción libre pero (1) `Tint` constante (no se actualiza) y (2)
acumula `Qcool/Qheat` en vez de `Qin/Qout`. Devuelve
`(Ti_series(=setpoint), Tso, Tsi, [Th], Tfield, days, Qcool, Qheat)`.

**API (`eh2d.System2D`), espejo de `System`:**
- `solveAC()` → corre el día con AC; devuelve `Ti` (Series **constante** = setpoint), guarda
  `cooling_energy=Qcool`, `heating_energy=Qheat`, `energy_transfer=None`. Caché separada del
  `solve()` libre (la firma incluye el modo AC, como `__last_solve` en el 1D).
- Setpoint: default `Tn.mean()` (1D); opción de atributo para fijarlo.
- `info()` ya imprime placeholders de `cooling/heating energy`.

**Pruebas (`tests/test_eh2d_ac.py`):**
1. **Reduce al 1D:** muro/techo de capa homogénea con `solveAC` ≈ `System.solveAC` 1D
   (`atol` temperatura, `rtol ~1–2%` energía).
2. **Sanidad:** `Qcool≥0`, `Qheat≥0`; clima cálido → `Qcool>0`; `Ti` constante = setpoint.
3. **Periodicidad:** converge antes de `max_days`.
4. **AIRE vs RELLENA:** el hueco cambia la carga (mayor decremento → menor `Qcool`).

**Notas:** acoplable con `dt`/serial-paralelo (extras). El setpoint fijo es el modelo del 1D;
una banda de confort (calienta < Tlow, enfría > Thigh, flota en medio) sería una extensión
posterior, no parte de este espejo.

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
