# HALF-NODES: nodos en superficie con medios volúmenes de control

Documento de diseño. Cómo pasar del esquema actual (malla centrada en celdas,
volúmenes completos en todas partes) al esquema clásico de Patankar con nodos
sobre las superficies y medios volúmenes de control en las fronteras en contacto
con aire, en el 1D y en el 2D. No está implementado; este documento registra el
diseño, el impacto esperado y el orden de trabajo.

## 1. Situación actual (verificada en el código)

- 1D (`ehtools.set_k_rhoc`): `dx = L/Nx`, `Nx` celdas completas que cubren
  `[0, L]`; el nodo está en el centro de cada celda. La masa térmica es
  `ρc·dx/dt` para todas las celdas (`prepare_static_coefficients`), incluidas
  las dos de frontera. La película entra como conductancia directa sobre la
  celda de frontera: `a[0] = masa + ho + G`, `a[Nx-1] = masa + G + hi`.
- 2D (`eh2d.compute_mesh*`, `ehtools2d`): `dx = X/nx`, `dy = Y/ny`, celdas
  completas; `a_P = ρc·dx·dy/Δt` uniforme; las fronteras con aire (`j = 0`,
  `j = ny-1`) reciben `h·dx` sumado al nodo de la celda de frontera.

Consecuencias del esquema actual:

1. El "nodo de superficie" está a media celda del interior del sólido: la
   resistencia de conducción `(Δ/2)/k` entre la superficie real y el nodo no
   está en serie con `1/h`. Error de primer orden, ~0.3 % con las mallas por
   defecto (verificado contra la solución analítica Dirichlet: 0.58 % a
   `Nx = 200`, 0.13 % a `Nx = 800`).
2. `Tso`/`Tsi` reportadas no son temperaturas de superficie exactas sino del
   nodo a media celda.
3. La convención de reporte `/(nx-1)` heredada del C ("half nodes") presupone
   una malla que ya no existe; en la malla actual introduce el sesgo
   `+T/(nx-1)` documentado durante la validación.

El esquema de medios volúmenes elimina 1 y 2, y da sentido geométrico a 3
(aunque el promedio correcto es el ponderado, ver §4).

## 2. Implementación 1D

### 2.1 Malla

- Nodos en `x_i = i·Δx`, `i = 0 … Nx-1`, con `Δx = L/(Nx-1)`. Los nodos 0 y
  `Nx-1` quedan exactamente sobre las superficies exterior e interior.
- Volumen de control del nodo `i`: `[x_i - Δx/2, x_i + Δx/2] ∩ [0, L]`, es
  decir `Δx/2` en las dos fronteras y `Δx` en el interior. La suma de
  volúmenes es `(Nx-1)·Δx = L`, exacta.

### 2.2 Cambios en `ehtools.set_k_rhoc`

- `dx = L_total / (nx - 1)` en lugar de `L_total / nx`.
- El mapeo material→celda ya trabaja con extremos arbitrarios (integra ρc y
  `1/k` sobre spans): basta con recorrer los volúmenes recortados
  `[max(0, x_i - dx/2), min(L, x_i + dx/2)]`. La conservación de masa
  `Σ ρc_j·L_j` se mantiene exacta por el mismo promedio ponderado por espesor.
- `Gf[f]` pasa a ser la conductancia entre nodos consecutivos, integrando
  `∫ dx'/k(x')` sobre `[x_f, x_{f+1}]` (spans nodo-a-nodo; misma maquinaria,
  distintos límites).

### 2.3 Cambios en `prepare_static_coefficients`

Único cambio real: masa de los nodos extremos a la mitad.

```python
mass_coeff = rhoc_array * (dx / dt)
mass_coeff[0] *= 0.5
mass_coeff[nx - 1] *= 0.5
```

Las ecuaciones de frontera conservan su forma (`a[0] = masa₀ + ho + G₀`), pero
ahora `ho` conecta el aire con un nodo que ESTÁ en la superficie: la película y
la conducción quedan en serie de manera natural y `T[0]`, `T[Nx-1]` son las
temperaturas de superficie exactas.

### 2.4 Lo que no cambia

- `calculate_coefficients`, `solve_PQ`, `solve_PQ_AC`: leen `mass_coeff` y las
  diagonales; no saben de la malla.
- El acoplamiento con el aire interior (`capacitance_factor`, criterio de
  `hi_flow` con `T[Nx-1]`) no cambia de forma; de hecho mejora, porque
  `T[Nx-1]` pasa a ser exactamente `Tsi`.
- La contabilidad de energías (`hi·(T[Nx-1] - Ti)·dt`) tampoco cambia de forma.

## 3. Implementación 2D

### 3.1 Malla y pesos

- `Δx = X/(nx-1)`, `Δy = Y/(ny-1)`; nodos sobre los cuatro bordes del dominio
  (superficies exterior/interior en `j = 0, ny-1`; planos adiabáticos en
  `i = 0, nx-1`).
- Pesos de volumen por dirección: `w_i = ½` si `i ∈ {0, nx-1}`, si no `1`;
  igual `w_j`. Volumen del nodo: `(w_i Δx)(w_j Δy)` — medios volúmenes en los
  bordes, cuartos en las cuatro esquinas.

### 3.2 Cambios en el ensamblaje (`ehtools2d`)

En `_step_inner` / `_step_hueca` / `_step_slab` (y en `calculate_coefficients_2d`):

- Término de masa: `apo = ρc · (w_i Δx)(w_j Δy) / Δt`.
- Áreas de cara: la conductancia N/S usa el ancho de cara `w_i·Δx`; la E/W usa
  `w_j·Δy` (las medias armónicas de `k` no cambian de forma).
- Películas de frontera: en `j = 0`, `a_P += ho·(w_i Δx)` y
  `d += ho·(w_i Δx)·Tsa`; en `j = ny-1`, lo mismo con `hi`. La suma de anchos
  `Σ w_i Δx = X` es exacta.
- Flujo al aire interior (los 6 drivers):
  `flux = Σ_i hi·(w_i Δx)·(T[i, ny-1] - Ti)` — sigue siendo la integral exacta.
- Lo más simple es precalcular un arreglo `wx[nx]` (y `wy[ny]`) y pasarlo a los
  kernels, en lugar de ramificar por índice dentro de los bucles.

### 3.3 La cavidad: mantenerla como está (primera fase)

Los muros de la cavidad son fronteras internas cuya posición depende de la
geometría del usuario; no se puede garantizar que caigan sobre nodos con malla
uniforme. Propuesta pragmática: en la primera fase los medios volúmenes se
introducen SOLO en las cuatro fronteras externas del dominio; el tratamiento de
la cavidad (caras ajustadas a líneas de malla por redondeo, como hoy) no se
toca. Extender nodos-en-cara a las paredes de cavidad es una segunda fase
independiente y de mucho mayor alcance (cambia el carve-out, los índices
`i1,i2,j1,j2`/`cav_*` y los tipos de nodo NT).

### 3.4 Reporte de superficies y criterio de `hi_flow`

Con nodos en superficie, el promedio de superficie correcto es el ponderado:

```
Tsi = ( Σ_i w_i · T[i, ny-1] ) / (nx - 1)        # Σ w_i = nx - 1
```

- Reemplazar la convención actual `ΣT/(nx-1)` por este promedio ponderado en
  `Tso_series`/`Tsi_series` (nota: incluso en el C original, la suma plana
  entre `nx-1` era una aproximación del ponderado; el sesgo era
  `(T₀+T_last)/(2(nx-1))`).
- El criterio de `hi_flow` usa el mismo promedio ponderado.

## 4. Impacto esperado y compatibilidad

- Magnitud: corrimientos de primer orden, del orden del error de frontera
  actual (≲ 1 % en energías con mallas por defecto). `Tso`/`Tsi` cambian algo
  más (dejan de estar a media celda y pierden el sesgo de reporte).
- Golden tests del C (`legacy=True`): no tocar esos caminos; siguen ejecutando
  el esquema actual para la regresión de fidelidad.
- Pre-cómputos de la documentación (`docs/run_examples.py`): regenerar los
  tres casos tras el cambio.
- `numerics.qmd`: actualizar la ecuación discreta (`a_P` con volumen `w_iΔx·w_jΔy`)
  y la nota de mapeo capa→malla; `api.qmd`: la convención de muestreo de
  `Tso`/`Tsi`.

## 5. Criterios de aceptación

1. **Dirichlet analítico (1D)**: el script de validación por matrices de
   transferencia (AC, `ho = hi = 1e6`, `absortance = 0`) debe converger al
   valor analítico con error menor que el actual a igual `Nx` (hoy: 0.58 % a
   `Nx = 200`). Con nodos en superficie el límite Dirichlet es exacto en la
   frontera, así que el error restante debe ser solo de discretización interior.
2. **Conservación de masa térmica (1D)**: `Σ vol_i·ρc_i = Σ ρc_j·L_j` exacta a
   cualquier `Nx` (test existente; debe seguir pasando con los volúmenes
   recortados).
3. **Reducción 2D → 1D**: el `HollowBlock` relleno homogéneo debe reproducir el
   1D (test existente `test_reduces_to_1d`); ambos esquemas deben cambiar de
   forma consistente para que la equivalencia se mantenga.
4. **Cierre de energía**: `energy_imbalance ≈ 0` en régimen periódico (test
   existente).
5. **Golden C intactos**: los caminos `legacy` no cambian ni un bit.

## 6. Orden de trabajo sugerido

1. 1D completo (`set_k_rhoc` + `prepare_static_coefficients`) + aceptación 1 y 2.
2. 2D fase 1 (pesos `wx`/`wy` en masa, caras, películas y flujos; reporte
   ponderado y criterio `hi_flow`) + aceptación 3 y 4.
3. Regenerar pre-cómputos de docs y actualizar `numerics.qmd`/`api.qmd`.
4. (Opcional, fase 2) nodos sobre las paredes de cavidad.
