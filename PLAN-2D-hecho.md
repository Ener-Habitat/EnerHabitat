# PLAN-2D — Bitácora de lo hecho (Fases 0–8a)

Detalle de las fases completadas del port 2D (objetivo, entregables, archivos y
resultados de prueba contra el C). El plan vivo (diseño + pendiente) está en
[`PLAN-2D.md`](PLAN-2D.md).

## Fases

Cada fase es un PR pequeño con su prueba. Una fase no se da por buena hasta que su prueba pasa.

### Fase 0 — Golden master del C legacy ✅ HECHA
**Objetivo:** compilar y correr el C en esta máquina, congelar referencias.
- ✅ `legacy_eh/2dTfree/standalone/` con el C parchado: stubs vacíos de
  `database_begin/database_end` (eran las únicas consumidoras de `libpq/PQ*`), eliminado el
  `system("java … Mailparametros")` y el `#include <libpq-fe.h>`. `Makefile` nuevo sin `-lpq`,
  con `readline/ncurses` de Homebrew y `-std=gnu89` (el C usa retornos *implicit-int* de K&R
  que clang moderno rechaza). `dt=1`, `nx=ny=160` del `.inp` intactos.
- ✅ Flag de *debug dump* (`make DUMP=1`, `-DDUMP`) que escribe `dat/dump_meta.dat` (nx,ny,X,Y,
  dx,dy,i1,j1,i2,j2,tipo), `dump_NT.dat`, `dump_k.dat`, `dump_rhoc.dat` (formato `i j valor`),
  una vez tras armar la geometría. (`a,b,c,d`/`T` en un paso fijo quedan para cuando los pida
  la Fase 2/3.)
- ✅ Corrido el caso del `.inp` (bovedilla rellena, `tipo 2`). Series (`gbv_5_1.csv` **con
  columna de tiempo `t[h]`**), índices (`indice_gbv_5_1.csv`) y dumps congelados en
  `tests/golden/2d/` (ver su `README.md`).
- **Resultado:** binario compila y corre. **Baseline de tiempo: 722.55 s/día**
  (`nx=ny=160`, `dt=1 s`, macOS arm64, `-O2`).

**Cambios de alcance acordados durante la Fase 0:**
1. **Columna de tiempo en el CSV.** Se activó (estaba en el `printf` comentado) la escritura
   de `t/3600.` como **primera columna**, en **horas decimales** (1.0 = 1 h, 1.5 = 1 h 30 min).
   Encabezado nuevo: `t Is Tsa Ta T_paredext Tparedint Tint T_n DeltaT_n` / `[h] [W/m2] [oC]…`.
   Antes esta columna se pegaba desde otro archivo.
2. **Corrida de un día (configurable).** Se añadió un tope `MAXDAYS` al lazo de convergencia
   día‑a‑día: `make DAYS=N` corre N días (`DAYS=1` = un día y para); sin `DAYS` (o `0`)
   converge a `error<1e-5` como el original. El golden de Fase 0 es de **un día** desde la
   condición inicial uniforme, **no** el régimen convergido → **no** coincide con el
   `dat/gbv_5_1.csv` committeado (ése es convergido, sin columna de tiempo). C y Python
   correrán el mismo único día para validar. El golden convergido completo queda disponible
   con `make` sin `DAYS` si se necesita como referencia de régimen.

### Fase 1 — Geometría/topología (sin física) ✅ HECHA
**Objetivo:** `eh2d.py`: construir `NT`, `k`, `rhoc`, `dx,dy,X,Y,i1,j1,i2,j2` para bovedilla
rellena (`tipo 2`).
- ✅ `src/enerhabitat/eh2d.py`: enum `Bovedilla` (+`TIPO_C`), `compute_mesh` (cálculo de
  malla del `main`: `X,Y,dx,dy,i1,j1,i2,j2` con el mismo truncamiento entero del C),
  `draw_rellena` (port de `draw_viguetabovedilla2rellena`), `set_krhoc_rellena` (port de
  `set_krhocrelleno`), dataclass `Section2D.build()`, e inspectores `print_node_scheme` /
  `plot_node_scheme`. Matrices con forma `(nx,ny)` indexadas `[i,j]` (igual que el dump).
- ✅ **Prueba** `tests/test_eh2d_geometry.py`: lee el `.inp` + dumps de la Fase 0 y compara
  `NT` (igualdad exacta), `k`/`rhoc` (`rtol 1e-12`) y `X,Y,dx,dy,i1,j1,i2,j2`. **Pasa** en la
  malla `nx=ny=160` del golden. Verificado además el esquema visual en malla chica.
  *(Nota: `pytest` no está instalado en `.venv`; la prueba es auto-ejecutable como script:
  `.venv/bin/python tests/test_eh2d_geometry.py`.)*
- *Pendiente menor:* el cálculo de **índices** del `main` (Qin, decremento, retardo, …) se
  porta junto con el lazo temporal en la Fase 4, donde se necesitan.

### Fase 2 — Ensamble de coeficientes ✅ HECHA
**Objetivo:** `ehtools2d.calculate_coefficients_2d` → `a,b,c,d` para un `T` dado (bovedilla
rellena: solo casos `NT` 1–8 y 13).
- ✅ `src/enerhabitat/ehtools2d.py`: `calculate_coefficients_2d(NT,k,rhoc,To,T,Tsa,Tint,ho,
  hi,dt,dx,dy)`. Port vectorizado (numpy) de `calculate_coefficients` para `NT ⊆ {1-8,13}`;
  media armónica en caras, vecinos N/S diferidos a `d`, `a=aP`, `b=aE`, `c=aW`. Lanza
  `NotImplementedError` si aparecen nodos de hueco (0,9-12). (JIT numba → Fase 5.)
- ✅ **Dump del C** (añadido al standalone): `make DUMPCOEF=1` llena `T,To` con un campo
  determinista (rampa + ruido entero, reproducible exacto en Python), `Tsa=30,Tint=24`, llama
  `calculate_coefficients` una vez y vuelca `dump_coef_{a,b,c,d,T,To}.dat` + `dump_coef_meta`
  (Tsa,Tint,ho,hi,dt), luego sale (sin lazo temporal). Congelados en `tests/golden/2d/`.
- ✅ **Prueba** `tests/test_eh2d_coeffs.py`: compara `a,b,c,d` con el dump usando el mismo
  campo. **Pasa** con error relativo ~1e-16 (≪ `rtol 1e-10`). Muestra nodos representativos
  (vigueta vs bovedilla) lado a lado Python/C.

### Fase 3 — Solver de un paso (línea-TDMA + aire interior) ✅ HECHA
**Objetivo:** `ehtools2d.solve_PQ_2d`: barrido Gauss-Seidel por líneas con lazo interno
`|error|>1e-10` + actualización lumped de `Tint`.
- ✅ `src/enerhabitat/ehtools2d.py`: `solve_step_2d(...)` + helper `_tdma_rows` (TDMA en x
  vectorizado sobre todas las filas, mismo orden de operaciones del C). Lazo interno
  `do…while(|error|>1e-10)` con `error = Σ(T-Tn)/T/(nx·ny)` (con signo), coeficientes
  recalculados cada iteración con la `T` más reciente (vecinos en y diferidos), y
  actualización del aire interior `Tint = (Qh + (ρc·La·X/dt)·Ti)·dt/(ρc·La·X)`.
- ✅ **Dump del C** (`make DUMPSTEP=1`): un `solve_PQ` desde un campo determinista, vuelca
  `dump_step_{T0,T}.dat` + meta (Tsa,Tint_in/out,ho,**hi efectivo**,dt,La,X,ρair,cair,Qin,
  inner_iters). Congelado en `tests/golden/2d/`.
- ✅ **Hallazgo:** `convective_coefficients` **sobrescribe `hi`** según `beta` (muro `beta=90`
  → `hi=8.1`, ignora el `hi=6.6` del `.inp`); el solver usa 8.1. Se replica. `interchange(A,B)`
  copia `A←B`.
- ✅ **Prueba** `tests/test_eh2d_step.py`: **iteraciones internas py=32 = c=32**; campo `T`
  max|Δ|=**1.4e-14 °C**; `Tint` coincide exacto — todo ≪ `atol 1e-6`. Muestra el perfil de
  temperatura exterior→interior.

### Fase 4 — Integración temporal + convergencia día-a-día (bovedilla rellena) ✅ HECHA
**Objetivo:** lazo `t∈[0,86400]` y `while(error>1e-5)`, alimentado por la Tsa sinusoidal
del C (portada al harness de prueba). Calcular series e índices.
- ✅ **Frontera del C** portada a `tests/c_boundary.py`: `sol.h` (día juliano, declinación,
  orto/ocaso, duración del día) + `Ta_Tc_DtaT` + `time_evolution_Tsa`. Validada contra las
  columnas `Tsa,Ta,Is` del CSV golden (max|Δ| ≤ redondeo a 2 decimales; `Tc=29.21`,
  `DtaT=1.25` exactos).
- ✅ **Driver del día** `solve_day` (en `tests/test_eh2d_fullday.py`): replica el `main` para
  bovedilla rellena, **un día** desde IC uniforme (igual que el golden `DAYS=1`). Detalles
  finos replicados: `Tso` pre-solve `/nx`, `Tsi` post-solve `/(nx-1)`, `hi=8.1` (muro),
  acumuladores TPI/DDH/Qin, índices con `contador=86400/dt`. Guarda caché
  `golden/2d/py_day_5_1.csv` + `py_indice_5_1.csv`.
- ✅ **Prueba** `tests/test_eh2d_fullday.py` (~28 min de corrida, cacheada): serie
  `Tsa,Tso,Tsi,Tint` **max|Δ|=0.005 °C** (= redondeo del golden) ≪ `atol 0.1`; **10 índices**
  dentro del redondeo (decremento 0.43, retardo 6.1 h, TPI 76/68, DDH 34.5/14.7). Visual:
  tabla de índices y curva `Tsa(t)`/`Tint(t)` (amortiguamiento + retardo).
- *Nota:* solo se corre 1 día (golden `DAYS=1`); la convergencia día-a-día completa y el
  costo se atacan con JIT en Fase 5/7.

### Fase 5 — Integración al paquete (API 2D, producción) ✅ HECHA
**Objetivo:** API de alto nivel reutilizando EPW+pvlib.
- ✅ `eh2d.py`: clase **`System2D`** que reúsa la cadena clima/solar del paquete (la `Tsa(t)`
  sale de un `System` 1D interno → EPW+pvlib, al paso `config.dt`), arma la geometría por
  capas (+ bovedilla rellena opcional) y corre el motor del día. Exportada en `__init__`
  junto a `Section2D`, `Bovedilla`, `config2d`.
- ✅ **`config2d`** (en `config.py`): `nx, ny, tol_inner, tol_day, max_days` sin tocar el 1D.
- ✅ **JIT numba** (`ehtools2d.solve_day_2d` + `_step_inner`, `@njit(cache=True)`): motor del
  día completo con convergencia día-a-día, kernels compilados (≈ del orden de segundos/día tras
  compilar, vs ~28 min en numpy puro).
- ✅ **Correcciones conscientes** (respaldadas por la regresión de Fase 4, documentadas en
  `ehtools2d`): (1) actualización del aire interior con **un solo `dt`** (el C tenía un `dt²`
  latente que solo coincide a `dt=1`) → el 2D reduce EXACTO al 1D a cualquier `dt`;
  (2) promedios de superficie `Tso`/`Tsi` a `/(nx-1)`.
- ✅ **Prueba** `tests/test_eh2d_package.py` sobre el EPW de `tests/`:
  **reduce al 1D** (capa homogénea, `max|Δ|=3.1e-3 °C` ≪ 0.1) · **periodicidad** (converge en
  3 días) · **balance de energía** (`Qin=Qout`, desbalance 0.00 %). Visual: curvas `Ti(t)`
  1D y 2D superpuestas.

### Fase 6 — Bovedilla con cámara de aire (`tipo 1`) ✅ HECHA
**Objetivo:** radiación entre paredes (factores de vista), `hh` (Nusselt muro/techo),
nodo de aire del hueco `Thueco`, y casos `NT` 0, 9–12.
- ✅ **Geometría** (`eh2d`): `draw_hueca`, `set_krhoc_hueca`, `Bovedilla.AIRE` en
  `Section2D.build`. Golden `tipo 1` regenerado (`.inp` `tipo 1`, `cont=2`) en
  `tests/golden/2d/hueca/` (geometría + un paso + día).
- ✅ **Física del hueco** (`ehtools2d`): kernel njit `_step_hueca` con factores de vista
  (`_view_factors`), radiación Stefan-Boltzmann entre las 4 paredes, Nusselt
  `hh=0.4005·|ΔTud|^0.3033/e22^0.0901` (muro), nodo de aire `Thueco` (lumped, un solo `dt`),
  y casos `NT` 0/9-12. Wrappers `solve_step_hueca` (un paso) y `solve_day_hueca` (día, njit).
- ✅ **Prueba** `tests/test_eh2d_hueca.py`:
  - geometría exacta (0 diffs; 14124 aire, paredes 9/10=132, 11/12=107);
  - **un paso a precisión de máquina**: `hh` Δ=3e-16, `Thueco`/`Tint` exactos, campo `T`
    1.8e-14 °C, iters 32=32;
  - **día completo** (1 día, casa con golden `DAYS=1`): serie `Tsa/Tso/Tsi/Tint` max|Δ|≤0.013 °C
    e índices dentro del redondeo (decremento 0.65, retardo 4.0 h, TPI 63.9/49.3, DDH 51.8/23.5).
- *Notas:* sólo muro (`beta=90`); el techo (`beta=0`, Rayleigh) queda pendiente. `System2D`
  aún enruta al motor de relleno; conectar `Bovedilla.AIRE` a `solve_day_hueca` en producción
  es un seguimiento menor.

### Fase 7 — Desempeño y paralelización ✅ HECHA
**Objetivo:** medir y analizar viabilidad.
- ✅ **Hallazgo clave:** el método del C **ya es Jacobi por líneas** (en cada iteración interna
  se ensambla con un único snapshot de `T` y se resuelven todas las filas a un buffer antes de
  `T←Tnew`). Las filas son **independientes** → paralelizar el barrido sobre `j` con `prange`
  ejecuta el MISMO algoritmo, **sin costo de convergencia**. *Red-black es innecesario* (no hay
  acoplamiento Gauss-Seidel que desenredar).
- ✅ **Variante paralela** (`ehtools2d`): `_step_inner_par` + `solve_day_2d_par`
  (`@njit(parallel=True)`, `prange` en ensamble/TDMA/reducción; buffers `P,Q` locales por hilo).
  Portable: threading layer interno de numba (workqueue), **sin libs externas**.
- ✅ **Benchmark** `tests/test_eh2d_perf.py` (esta máquina, 18 núcleos):
  - **correctitud**: paralelo reproduce serial **al bit** (Δ=0.0);
  - **speedup por hilos modesto**: 1.27× @ 4 hilos, no escala (filas cortas, granularidad fina);
  - **escalado de malla** (serial): 0.1/0.5/6.0/19.6 ms/paso para 40²/80²/160²/240²;
  - **barrido de `dt`** (esquema implícito): subir `dt` da mejora **sublineal** porque el lazo
    interno necesita más iteraciones por paso — `dt=10` → ~1.6× con Δ<0.001 °C (punto dulce);
    `dt=60` Δ≈0.005 °C; `dt=300` Δ≈0.02 °C.
- **Recomendación final:** la paralelización con `numba prange` es portable y correcta pero
  rinde poco (~1.3×) por el TDMA secuencial en filas cortas; la palanca práctica mayor es
  **`dt=10`** (~1.6×, error despreciable). Combinadas ~2×. Ganancias mayores exigirían cambio
  algorítmico (line-Gauss-Seidel real para reducir iteraciones internas) o GPU (no portable).

---

## Fase 8a — Bloque hueco de concreto (muros): entregado

- `ehtools2d.solve_day_hueca_prod`: versión de **producción** del día con cámara de aire
  (aire interior con un solo `dt`, superficies `/(nx-1)`, + `Qin/Qout`); no toca la fiel.
- `eh2d.HollowBlock(material, emissivity, geometry)`: elemento 2D de muro (bovedilla `AIRE`,
  `required_tilt=90`, geometría amistosa `web/block_width/cover_*` + alias `a*/e*`). Exportado.
- `eh2d.System2D` **reescrita** (patrón de `System`): `Tsa()` por composición de un `System`
  1D, `_build_section()` mapea `layers`(+elemento) → `Section2D`, `solve()` enruta por tipo
  (AIRE→`solve_day_hueca_prod`), devuelve `pandas.Series` `Ti` alineada a `Tsa()`, guarda
  `energy_transfer/Qout/days/solve_dataframe`. Validación elemento↔`tilt`. Se borró el
  `_build_fields` simplificado.
- Fixture `tests/materials_2d.ini` (Concreto, Mortero, Yeso, EPS).
- **Prueba** `tests/test_eh2d_hollowblock.py` (5/5 pasan, malla chica): metodología (flujo
  idéntico al 1D), periodicidad (4 días), balance de energía (`Qin=Qout`, 0.00 %), guardas de
  orientación (`tilt≠90` falla) y de "un elemento". Demo: factor de decremento 0.47.

## Inspector a escala de la asignación de materiales: entregado

Implementado en `eh2d.py` y `System2D`: `section()`, `preview(field/panels/backend/save)`
(matplotlib a escala — ejes en mm, `aspect="equal"`, fronteras de capa, contorno del bloque,
leyenda; ASCII a escala de respaldo), `section_report()` (tabla de NT + materiales con `k,ρc`
y rango en `y`). matplotlib como extra opcional `enerhabitat[viz]`. Prueba
`tests/test_eh2d_inspect.py` (5/5).
