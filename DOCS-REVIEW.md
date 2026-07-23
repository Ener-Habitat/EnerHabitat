# Revisión crítica y propuesta de mejora de la documentación de EnerHabitat

- **Fecha de la revisión:** 22 de julio de 2026
- **Versión revisada:** EnerHabitat 0.2.1
- **Commit local revisado:** `66485d7`
- **Alcance principal:** `docs/`, `README.md` —que también se publica en PyPI—,
  `CITATION.cff`, implementación en `src/enerhabitat/` y pruebas en `tests/`. El
  manuscrito de `softwareX/EnerHabitat/` se consultó únicamente para detectar
  afirmaciones incompatibles entre el artículo y la documentación.

## Dictamen

La documentación tiene una estructura clara y ya declara varias limitaciones
importantes del modelo. Sin embargo, **todavía no debería presentarse como
documentación científica definitiva ni usarse para sostener, sin reservas, que la
versión 0.2.1 está validada**. Hay afirmaciones que contradicen la implementación,
una atribución normativa incorrecta, una formulación radiativa que no corresponde
al supuesto de superficies grises y criterios numéricos que no garantizan la
convergencia que el texto afirma.

La recomendación es resolver primero los asuntos P0 de este documento y después
reescribir las páginas. Varios de ellos requieren una decisión o una corrección en
el código: no se pueden arreglar solamente suavizando la narrativa.

Esta revisión no constituye una certificación formal del modelo ni sustituye una
revisión por pares. Es una auditoría de trazabilidad entre documentación,
ecuaciones, código, pruebas y fuentes primarias.

## Escala de prioridad

- **P0 — bloqueante:** puede cambiar resultados, su interpretación científica o la
  validez de una afirmación central.
- **P1 — mayor:** no necesariamente invalida el cálculo, pero impide reproducirlo o
  delimitar correctamente su dominio de aplicación.
- **P2 — moderada:** error de API, consistencia, terminología o reproducibilidad que
  debe corregirse antes de publicar.
- **P3 — editorial:** estilo, navegación o mantenimiento.

## Fortalezas que conviene conservar

- La documentación distingue los modelos 1D y 2D y organiza por separado física,
  método numérico, uso y API.
- El alcance de «un solo componente opaco» y la exclusión de ventanas, ventilación,
  infiltración y ganancias internas ya aparecen cerca del inicio.
- Se indican el sentido exterior → interior, las unidades de propiedades y la
  normalización de energías por área.
- Se reconoce que la versión Python difiere de la herramienta publicada en 2016;
  esa tabla es una buena base, aunque obliga a separar también sus validaciones.
- Las pruebas de regresión contra el código C, las pruebas de reducción 2D → 1D y
  los ejemplos versionados son activos valiosos. Deben describirse con la categoría
  de evidencia que realmente aportan.

## Bloqueantes científicos (P0)

**Todos atendidos** (jul 2026): P0-01, P0-03, P0-04 y P0-05 corregidos en este
repo; P0-02 resuelto con radiosidades (ver `RADIOSIDAD.md`); P0-08 cerrado por
decisión editorial. **P0-06 y P0-07** (matriz de evidencia, estudios de malla /
paso / tolerancias, y validación contra EnergyPlus y el experimento de Borbón)
se atienden en el **repositorio de validación**, aparte de este.

## Homogeneización de las páginas de teoría (plan, jul 2026)

Análisis comparativo de `model-1d.qmd` y `model-2d.qmd`: no comparten
estructura. Desajustes detectados:

1. El 2D **no tiene lista de supuestos** (sección invariante fuera del plano,
   aire de cavidad bien mezclado/no participante, un `h_c` uniforme por
   cavidad evaluado con temperaturas medias, superficies grises a temperatura
   media, lados adiabáticos) — el 1D sí la tiene.
2. `model-2d:16` «By periodicity, no heat flows» — incorrecto: la
   periodicidad no impone flujo cero; lo adiabático vale porque los lados son
   **planos de simetría** (hallazgo de la auditoría 2D).
3. Modos de solución invisibles en 2D: falta declarar el flujo **promediado
   en anchura** que alimenta al nodo de aire y que en AC **las cavidades
   siguen flotando**.
4. `### Solid fill (Fill.SOLID)` anidado bajo `## Cavity physics (Fill.AIR)`.
5. El 2D no tiene tabla de *Outputs and units* (el 1D sí) y omite sus salidas
   extra (`Tso/Tsi/Thueco/Tfield/Qout/solve_dataframe`).
6. Títulos no paralelos («physical problem» vs «non-homogeneous systems»).
7. El mapeo de ejes 1D↔2D ($x$ del 1D = $y$ del 2D) no está declarado.

**Esqueleto común acordado**: Intro → *Domain and assumptions* (2D: ejes +
mapeo + supuestos como deltas del 1D) → *Governing equations* → *Boundary
conditions* (1D: sun–air + indoor/modes; 2D: bloque compacto con puntero) →
[2D: *Solution modes* (deltas), *Cavity physics* (con Solid fill como
hermano), *Geometries*] → [1D: *The average day*; 2D: puntero] → *Outputs and
units* (tabla en ambos; 2D añade Constraints) → *References*. Verificar
anclas (`model-*.qmd#...`) desde numerics/usage/api al renombrar secciones.

## Prioridades de la narrativa restante

**N1 — corrección científica del sitio (antes del envío):**

- [x] Homogeneizar `model-1d`/`model-2d` según el plan de arriba (incluye
  corregir simetría vs periodicidad, cavidades en AC, supuestos 2D). *(Hecho
  jul 2026.)*
- [x] Corregir los **errores factuales** del día promedio: `Ib` es **DNI**
  (normal directa), no «beam horizontal»; declarar modelo de transposición y
  albedo de pvlib (isotrópico, 0.25); semántica real de `day`/`year` en TMY;
  firma real de `meanDay()`. *(Hecho jul 2026; también declarada la
  acumulación horaria Wh/m² del EPW tratada como potencia media.)*
- [x] Tabla única de supuestos y consecuencias, enlazada desde la portada
  (se alimenta de la homogeneización). *(Hecho jul 2026:
  `docs/assumptions.qmd`, en el menú Theory, enlazada desde portada y ambas
  páginas de modelo.)*
- [x] Bibliografía: NOM con edición/apéndice, DOIs faltantes, procedencia del
  3.9 K (Mackey & Wright 1944 + ASHRAE), especificación EPW citada. *(Hecho
  jul 2026, salvo dos pendientes: nota de derivación de la correlación de
  Xamán — se confirmó algebraicamente que es la reducción dimensional de un
  ajuste `Nu = C·Ra^0.3033`, falta confirmar del autor cuál ajuste/regímen —
  y la decisión sobre la errata de `ν` (ver Hollands, arriba).)*
- [x] **Regenerar** los CSV 2D precalculados y cifras de ejemplos. *(Hecho
  jul 2026, con toda la física final: radiosidad + mapeo 1D + ν + C_w de
  Xamán. Malla 80×160: muro libre −0.20 %, muro AC Qcool −2.26 % / Qheat
  −3.66 %, techo libre −0.61 % vs publicados; usage.qmd re-renderizado.
  Pendiente del lado del manuscrito: cotejar las cifras 2D que cite.)*

**N2 — exactitud de afirmaciones y derivaciones:**

- [ ] Reescritura de afirmaciones centrales (Portada/README y texto
  provisional de validación mientras el repo de validación produce la matriz).
- [ ] Energías discretas: definir/exponer `Q_in`/`Q_out`/`Q_net` y la
  convención temporal de las salidas 1D (la del 2D ya está declarada en api).
- [ ] Derivación de la frontera discreta (media celda, `G_∞P`) con figura de
  malla; decidir si se corrige el esquema o solo se documenta el actual.
- [ ] Auditoría 2D restante: promedio superficial `/(nx-1)` de `Tso`/`Tsi` y
  dimensiones efectivas en `section_report()`.
- [ ] Método numérico: «TDMA exactly» → «salvo redondeo»; retirar «buying
  almost no speed» o respaldarlo.
- [ ] Temperatura neutral: precisar que `T_n` usa la media del día sintético;
  verificar la fuente de la tabla `DeltaTn` en Morillón 2004.

**N3 — editorial y mantenimiento:**

- [ ] Manifiestos de reproducibilidad en los ejemplos (versión, commit,
  hashes, malla, tolerancias, `converged`/`day_error`/`energy_imbalance`).
- [ ] Retirar o respaldar cifras de rendimiento («10–20 min», «≈6×», «≈10×»);
  procedencia y licencia del EPW incluido.
- [ ] PyPI: render de ecuaciones del README, reducir duplicación
  README↔Quarto, tono de «we love uv».
- [ ] Arquitectura documental (separar *Theory* de *Evidence*) — opcional
  antes del artículo; natural cuando llegue la matriz de validación.

## Auditoría de las ecuaciones 1D

### Ecuación de conducción

La ecuación

$$
\rho_jc_j\frac{\partial T_j}{\partial t}
=k_j\frac{\partial^2T_j}{\partial x^2}
$$

es correcta para una capa homogénea, isotrópica y con propiedades constantes, sin
generación volumétrica. La documentación debe agregar explícitamente: contacto
térmico perfecto, resistencia de contacto nula, ausencia de cambio de fase,
transporte de humedad y fuentes internas.

Para mantener una notación única, usar `c` o `c_p` en todo el proyecto; actualmente
`docs/` usa `c` y `README.md` usa `c_p`.

### Interfaz entre capas

`docs/model-1d.qmd:56-60` usa la misma etiqueta de evaluación en ambos lados y no
muestra explícitamente la continuidad de temperatura. Una formulación rigurosa es

$$
T_j(x_j^-,t)=T_{j+1}(x_j^+,t),
$$

$$
-k_j\left.\frac{\partial T_j}{\partial x}\right|_{x_j^-}
=-k_{j+1}\left.\frac{\partial T_{j+1}}{\partial x}\right|_{x_j^+}.
$$

Si se admite resistencia de contacto en el futuro, esta primera igualdad deberá
reemplazarse por el salto correspondiente.

### Condiciones exteriores e interiores

Los signos de las ecuaciones continuas son coherentes si `x` crece del exterior al
interior y se declara ese convenio de flujo. Falta, no obstante, derivar cómo se
aplican al volumen de frontera. El código añade `h_o` o `h_i` directamente al nodo
extremo y asigna a ese nodo una capacidad de volumen completo `ρcΔx`. Si los nodos
son centros de volúmenes, la conductancia superficie–centro debería considerar la
media celda:

$$
G_{\infty P}=\frac{A}{1/h+\delta_P/k_P}.
$$

Si son nodos sobre la superficie, la capacidad y las distancias deben reflejar
volúmenes de media anchura. La documentación actual no permite saber cuál es la
interpretación pretendida y afirma que los coeficientes entran «exactamente» como
las condiciones de frontera. Se requiere una figura de malla y una derivación,
seguida de una prueba de orden de convergencia.

### Temperatura sol–aire

La expresión

$$
T_{sa}=T_a+\frac{\alpha_s I_s}{h_o}-RF
$$

es una condición equivalente empírica útil, pero `RF` no es un cálculo explícito de
intercambio con cielo y entorno. *(Resuelto jul 2026: `RF` ahora varía linealmente
3.9 → 0 entre 0° y 90°, la regla de la herramienta de 2016; la docs lo declara como
corrección empírica.)* Pendiente:

- citar capítulo, edición, ecuación y condiciones de procedencia de 3.9 K;
- usar `α_s` para absorptancia solar y reservar `a` para coeficientes numéricos.

La API escribe `absortance`; debe identificarse como una grafía heredada de la API,
mientras el texto físico usa *absorptance*.

### Nodo interior libre y definición discreta de energía

La ecuación del nodo interior es dimensionalmente correcta **por cada metro cuadrado
de componente**: `L_a` equivale a un volumen de aire de `L_a m³` por `1 m²` de
superficie. Esto debe decirse expresamente; no es la profundidad de una habitación
completa salvo que su relación volumen/área coincida.

La integral positiva de `docs/model-1d.qmd:128-135` representa energía bruta que
entra al nodo, no transferencia neta diaria. El código también calcula internamente
la energía bruta saliente. En régimen perfectamente periódico, la neta debería
cerrar cerca de cero, pero el solver no usa ese cierre como criterio.

**Propuesta:** definir y exponer por separado

$$
Q_{in}=\int_0^P\max(q_i'',0)\,dt,
\quad
Q_{out}=\int_0^P\max(-q_i'',0)\,dt,
\quad
Q_{net}=Q_{in}-Q_{out},
$$

así como el instante temporal usado en la cuadratura. `energy_transfer` puede
mantenerse como alias histórico de `Q_in`, pero no debe llamarse neta. La igualdad
`Q_in = Q_out` debe mostrarse como verificación numérica con tolerancia, no como una
premisa.

### Temperatura neutral y banda de confort

El código calcula `T_n = 13.5 + 0.54 mean(T_a)` usando la temperatura **sintética del
día promedio**, no la media mensual cruda del EPW. La documentación debe precisar
esta diferencia.

La tabla escalonada de `DeltaTn` se atribuye a Morillón, pero la referencia actual no
incluye página, tabla o ecuación. El artículo
[Morillón et al. (2004)](https://doi.org/10.1016/j.solener.2003.11.008) debe
revisarse para confirmar que respalda exactamente esos diez intervalos. Si la tabla
procede de otra publicación o del software legado, hay que citar esa fuente
primaria. También debe evitarse presentar la banda como un criterio de confort
normativo universal.

## Construcción del día promedio y forzamiento solar

Esta parte necesita una especificación algorítmica, porque `meanDay()` no extrae un
«día meteorológico promedio» directo. La implementación actual hace lo siguiente:

1. Selecciona todas las filas del mes en el EPW.
2. Calcula las medias mensuales de los mínimos y máximos diarios de temperatura y la
   hora media del máximo diario.
3. Construye `T_a(t)` con el modelo de Chow–Levermore, colocando el mínimo al
   amanecer solar del día sintético elegido.
4. Promedia `Ig`, `Id` e `Ib` por hora del reloj dentro del mes.
5. Interpola linealmente esas 24 medias a 1 s; desde las 23:00 hasta 23:59 conserva
   el último valor mediante *forward fill*, sin cierre periódico 23:00 → 00:00.
6. Calcula posición solar para la fecha sintética y después submuestrea a 10 s.

### Nombres y unidades de irradiancia

La frase «global, beam and diffuse horizontal irradiance» en
`docs/model-1d.qmd:176-177` es incorrecta para `Ib`. El mapeo real del EPW es:

| Campo EnerHabitat | Campo EPW | Significado |
|---|---|---|
| `Ig` | Global Horizontal Radiation | GHI, irradiancia global horizontal, W/m² |
| `Ib` | Direct Normal Radiation | DNI, irradiancia directa normal, W/m² |
| `Id` | Diffuse Horizontal Radiation | DHI, irradiancia difusa horizontal, W/m² |
| `Is` | `poa_global` de pvlib | irradiancia total sobre el plano, W/m² |

La especificación EPW define los tres campos solares en Wh/m² acumulados durante
el intervalo precedente; por ejemplo, la hora 1 cubre 00:01–01:00
([diccionario EPW de EnergyPlus](https://bigladdersoftware.com/epx/docs/22-1/auxiliary-programs/energyplus-weather-file-epw-data-dictionary.html)).
En un intervalo de una hora, el valor numérico puede convertirse en potencia media
en W/m² dividiendo por una hora, pero no es una muestra instantánea. El código resta
uno a la hora, trata las medias como puntos horarios e interpola linealmente. Debe
justificarse esa reconstrucción, fijar si cada valor se coloca al inicio, centro o
fin del intervalo y comprobar que la integración de la señal reconstruida conserva
la energía solar del EPW. Esta comprobación también debe cubrir archivos con más de
un intervalo por hora.

El código llama a `pvlib.irradiance.get_total_irradiance()` sin especificar modelo
ni albedo. En pvlib 0.13.1 eso implica cielo isotrópico y albedo 0.25; además
`poa_global` suma componentes directa, difusa de cielo y reflejada por el suelo
([documentación oficial de pvlib](https://pvlib-python.readthedocs.io/en/v0.13.1/reference/generated/pvlib.irradiance.get_total_irradiance.html)).
Estos supuestos afectan el resultado y deben aparecer en teoría y API. Idealmente,
modelo de transposición y albedo deberían ser parámetros explícitos y formar parte
del manifiesto de cada simulación.

### Significado de `day` y `year`

En un TMY, `year=2025` **no selecciona observaciones de 2025**. La implementación
reemplaza el año de los registros para formar un índice compatible y usa la fecha
para geometría solar. `day` también elige la fecha solar/sintética; no selecciona
solo el día 15 del archivo. Esta semántica debe explicarse junto a la firma.

La firma real es

```python
meanDay(day="15", month="current", year="current")
```

no `meanDay(day=15, month, year)`, que ni siquiera es una firma Python válida con
parámetros requeridos después de uno predeterminado.

### Zona horaria y afirmación «cualquier lugar»

`Location` convierte el offset EPW con
`int(datos[8].split('.')[0])`. Esto trunca offsets fraccionarios. Un EPW con UTC
`+5.5`, `+5.75`, `+9.5`, etc. recibe una hora solar incorrecta. Mientras no se
corrija, «user-provided EPW file, any location» y «worldwide» son afirmaciones
falsas para una parte relevante de los archivos válidos.

**Propuesta:** conservar el offset decimal mediante una zona fija y añadir pruebas
para ±3.5, +5.5, +5.75 y −9.5. Documentar también que se usa tiempo local estándar
del EPW, no reglas históricas de horario de verano.

### Convención temporal de las salidas

La serie se etiqueta con el instante que alimenta `Tsa[s]`, pero `Ti[s]` se guarda
después de avanzar un paso. Esto introduce un desfase de un `Δt` en la interpretación
del sello temporal, aunque los índices coincidan. Debe elegirse y documentarse una
convención —estado al inicio o al final del intervalo— y aplicarse a `Tsa`, `Ti`,
`Tso`, `Tsi`, `Thueco` y flujos.

## Auditoría del modelo 2D

### Dominio sólido y fronteras laterales

La ecuación con `∇·(k∇T)` es adecuada para los sólidos con `k(x,y)` por tramos. Debe
aclararse que las celdas de aire no se resuelven con esa PDE: se sustituyen por un
nodo bien mezclado y leyes de intercambio en el perímetro.

«Por periodicidad no hay flujo» no es cierto en general. Una frontera periódica
impone igualdad de temperatura y continuidad de flujo entre lados opuestos; no
impone gradiente normal cero. El código aplica fronteras adiabáticas, apropiadas si
los lados elegidos son **planos de simetría** de la celda repetitiva. Reescribir
`docs/model-2d.qmd:14-17` y `46-54` en esos términos y mostrar por qué las geometrías
concretas sí tienen esa simetría.

### Nodo de aire de cavidad

La ecuación integral publicada representa bien el balance conceptual. Faltan cuatro
detalles del algoritmo:

- una temperatura uniforme por cavidad;
- un solo `h_c` uniforme en sus cuatro caras;
- `h_c` se obtiene con las temperaturas medias de las caras exterior/interior a
  través del espesor, no con cada temperatura local;
- `T_h` queda fija durante la iteración del sólido y se actualiza explícitamente al
  terminar el paso.

Por ello, la frase «the converged step is consistent with the cavity state» es
demasiado fuerte: la no linealidad del sólido se itera con el valor anterior de
`T_h`.

### Correlación de Xamán para muros

La fuente citada estudia cavidades rectangulares **altas**, con relaciones de
aspecto 20, 40 y 80 y números de Rayleigh de `10²` a `10⁸`
([Xamán et al., 2005](https://doi.org/10.1016/j.enbuild.2004.11.001)). El ejemplo
documentado de bloque tiene una cavidad de 0.16 × 0.08 m, muy lejos de esas
relaciones de aspecto según la convención usada.

Es necesario rastrear la procedencia exacta de

$$
h_c=0.4005|\Delta T|^{0.3033}/d^{0.0901},
$$

incluidas unidades, longitud característica, geometría, rango de `Ra` y tratamiento
laminar/turbulento. Si la expresión procede de una derivación intermedia de
Barrios et al. y no directamente de Xamán, se debe citar esa derivación. Hasta que
se demuestre aplicabilidad a bloques de baja relación de aspecto, la documentación
debe declararla como extrapolación.

**Resuelto (jul 2026):** rastreada a la **Ec. (11)** de Xamán 2005 (turbulenta,
A=20; el exponente de `d` = 3n−1 = −0.0901 prueba el linaje). La constante 0.4005
del C era una reducción no registrada ~0.61× de la fiel; **corregida**: producción
calcula `C_w = 0.0857·k·(gβ/να)^0.3033 ≈ 0.589` en runtime (propiedades a 300 K,
densidad configurable); el 0.4005 sobrevive solo en los golden de fidelidad al C.
Sensibilidad medida: +1 % en energías diarias del muro documentado (+4 % con la
variante pared–aire ×2, descartada por ahora). Extrapolación de relación de
aspecto declarada en `model-2d.qmd`. El caso Borbón de la campaña de validación
prueba constante y extrapolación directamente.

### Correlación de Hollands para techos

La forma publicada es la de una capa horizontal calentada desde abajo. La fuente
reporta datos de aire desde régimen subcrítico hasta `Ra = 4×10⁶`
([Hollands et al., 1975](https://doi.org/10.1016/0017-9310(75)90179-9)). Deben
documentarse y comprobarse:

- rango de `Ra` alcanzado por cada geometría permitida;
- definición y signo de `ΔT`;
- longitud característica `d`;
- condiciones de las paredes laterales y relación de aspecto;
- temperatura a la que se evalúan las propiedades.

El código usa constantes `k_air = 0.0262 W/(m·K)`,
`ν = 1.11×10⁻⁵ m²/s`, `β = 1/300 K⁻¹` y calcula `α` con la densidad y capacidad
configurables. Esa combinación debe tener fuente y temperatura de referencia; si el
usuario modifica propiedades puede dejar de ser termodinámicamente consistente.

**Hallazgo (jul 2026):** la fuente es Incropera (confirmado por el autor), y
`k_air`, `α` y `β` sí corresponden a ~300 K — pero **`ν = 1.11×10⁻⁵` no**:
Incropera a 300 K da `1.589×10⁻⁵`; el valor del código correspondía a aire a
~240 K (errata heredada del C). **Corregido (jul 2026):** ν ya no está
hardcodeada — se calcula `μ_Sutherland(300 K)/ρ_a` con la densidad
configurable (ν y α consistentes entre sí); test en
`tests/test_air_properties.py`. Solo afecta techos `Slab` AIRE. Pendiente:
verificación experimental con la campaña de validación.

### Geometría discreta y áreas

Las dimensiones se convierten a índices con truncamiento/redondeo (`int(...+0.5)`),
y `Δx = W/nx`, `Δy = e/ny`. Por tanto, la geometría calculada no coincide
exactamente con la solicitada. La documentación debe informar dimensiones efectivas,
error por segmento y número de celdas de cada material. `section_report()` es un
buen lugar para ello.

También debe revisarse el promedio de `Tso` y `Tsi`: varias rutas de producción
suman `nx` valores y dividen por `nx-1`, aunque la integración de flujo usa
`nx·Δx = W`. Si son volúmenes de ancho completo, el promedio consistente sería
`sum(T)/nx`; si existen medias celdas, se requieren pesos explícitos. Este punto
afecta las temperaturas superficiales guardadas, aunque no necesariamente las
energías que integran con `Δx/W`.

## Método numérico: redacción propuesta

La frase «fully implicit» debe limitarse a la discretización temporal de la
conducción sólida. El algoritmo completo es segregado:

1. conducción sólida con Euler hacia atrás;
2. solución TDMA exacta del sistema lineal 1D de ese paso;
3. en 2D, barridos línea por línea y punto fijo para términos en `y` y cavidad;
4. actualización explícita de `T_i` en modo libre;
5. actualización explícita de cada `T_h` tanto en modo libre como acondicionado;
6. repetición de días hasta el criterio parcial actual.

«TDMA solves exactly» debe decir «TDMA resuelve, salvo redondeo, el sistema
tridiagonal discretizado»; no significa solución exacta de la PDE. Del mismo modo,
la estabilidad incondicional de Euler implícito no implica exactitud temporal ni
estabilidad incondicional del esquema acoplado.

La media armónica

$$
k_f=\frac{2k_Lk_R}{k_L+k_R}
$$

solo es la forma directa para distancias iguales y contacto perfecto. Conviene dar
primero la forma por resistencias:

$$
G_f=\frac{A_f}
{\delta_L/k_L+R_c''+\delta_R/k_R},
$$

y mostrar cómo se reduce a la media armónica cuando
`δ_L = δ_R = Δx/2` y `R_c'' = 0`.

No debe afirmarse que bloquear `dt = 10 s` elimina la necesidad de un estudio
temporal. Un valor fijo puede ser estable y aun así introducir error de fase o
amplitud. Para una publicación, se necesita ejecutar el mismo caso con, por
ejemplo, 20, 10, 5, 2 y 1 s, o una secuencia equivalente, y reportar convergencia de
temperaturas, fase y energías. Si la API pública seguirá bloqueando `dt`, el estudio
puede hacerse mediante una ruta interna de verificación.

La afirmación de que aumentar `dt` «buying almost no speed» y que los barridos
crecen para compensar no está sustentada por resultados en la documentación; debe
eliminarse o acompañarse por un benchmark reproducible.

## Correcciones de API y narrativa (P1/P2)

**Todas atendidas** (jul 2026): tipos de retorno (Series), caché de `Tsa()`,
`config.file` con excepciones explícitas y rollback, `reset()` documentado,
diferencias `System2D`↔`System`, límite de 7 capas, validación dirigida de
rangos, `Fill.SOLID_SYMMETRIC` eliminado, salidas 2D documentadas,
diagnóstico de convergencia (vía P0-03/04), offsets UTC fraccionarios
(`pytz.FixedOffset`), «color»→absortancia, unidades SI unificadas y `RF`
lineal en la inclinación con capacidades 1D/2D separadas.

## Supuestos que faltan en una documentación científica

Incluir una tabla única de supuestos y consecuencias, enlazada desde la portada:

- materiales isotrópicos, homogéneos por región y con propiedades constantes;
- contacto perfecto y sin resistencia interfacial;
- sin humedad, difusión de vapor, condensación, lluvia, cambio de fase ni
  generación interna;
- coeficientes convectivos exteriores/interiores constantes;
- sin sombreado, entorno urbano, obstrucciones ni radiación de onda larga resuelta;
- nodo interior perfectamente mezclado y sin masa térmica de mobiliario u otras
  superficies;
- componente de área unitaria en 1D y sección invariante fuera del plano en 2D;
- laterales 2D como planos de simetría, no periodicidad general;
- aire de cavidad bien mezclado, no participante y con propiedades constantes;
- correlaciones de convección aplicadas fuera/dentro de su rango, según resulte de
  la revisión;
- setpoint constante e ideal, sin dinámica HVAC;
- día sintético mensual periódico, no una simulación anual ni una secuencia de días
  meteorológicos reales.

Cada salida debe enlazar los supuestos que más limitan su interpretación.

## Propuesta de reescritura de afirmaciones centrales

### Portada/README

Texto sugerido:

> EnerHabitat estima la transferencia sensible transitoria a través de un único
> componente opaco de la envolvente. Resuelve conducción 1D en sistemas
> multicapa y conducción 2D, con modelos reducidos de cavidad, en celdas
> constructivas específicas. El forzamiento es un día mensual sintético derivado
> de un EPW. Sus salidas son temperaturas y cargas térmicas ideales por unidad de
> área del componente; no representan por sí solas la demanda o el consumo de un
> edificio completo.

### Validación

Texto sugerido mientras no se complete la matriz propuesta:

> El modelo de la herramienta Ener-Habitat de 2016 fue comparado con EnergyPlus en
> 1D y con un experimento estacionario de bloque hueco en 2D. La versión Python
> 0.2.1 modifica entradas climáticas y parámetros numéricos. Sus pruebas locales
> verifican regresión respecto al código legado y consistencia entre
> implementaciones; estas pruebas no sustituyen una validación independiente de
> todas las configuraciones actuales. Los casos de validación reproducibles y su
> correspondencia exacta con cada versión se enumeran en [enlace versionado].

### Día promedio

Evitar «average EPW day» sin definición. Usar «día mensual sintético» y describir
en pseudocódigo los seis pasos enumerados antes. Añadir una figura que compare una
semana EPW, las medias horarias y la señal interpolada.

## Reproducibilidad de ejemplos y rendimiento

Los CSV 2D están precalculados y `freeze: auto` puede conservar resultados antiguos.
`summary.json` registra tiempo, días y energías, pero no la versión que los generó.
Cada ejecución publicada debe producir un manifiesto, por ejemplo:

```yaml
enerhabitat: 0.2.1
git_commit: 66485d7
python: 3.x.y
numba: x.y.z
pvlib: 0.13.1
numpy: x.y.z
pandas: x.y.z
epw_sha256: ...
materials_sha256: ...
date_parameters: {day: 15, month: 5, year: 2025}
mesh: {Nx: 200, nx: 80, ny: 160}
solver: {dt_s: 10, tol_inner: 1e-10, tol_day: 5e-4, max_days: 60}
convergence: {converged: true, day_error: ..., energy_imbalance: ...}
hardware: {cpu: ..., cores: ..., os: ...}
timing: {jit_warmup_excluded: true, repetitions: 5, statistic: median}
```

Las frases «10–20 minutes», «≈6× con 8 procesos» y «≈10× con 16» requieren una
tabla con caso, hardware, sistema operativo, versiones, compilación JIT fría/caliente,
número de repeticiones y dispersión. Si esos datos no están disponibles, retirar
las cifras y conservar solo una recomendación cualitativa.

El EPW incluido necesita fuente, fecha de descarga, periodo TMYx, licencia o
condiciones de redistribución y hash. Las propiedades «typical/illustrative» de
`docs/usage.qmd` deben incluir fuentes, temperatura/humedad de referencia e
incertidumbre, o mantenerse únicamente en un ejemplo claramente no destinado a
validación.

## Bibliografía y citación

`docs/references.bib` empieza con un comentario que dice que todavía se añadirán
NOM, ASHRAE, Duffie & Beckman y otras fuentes, aunque las NOM y ASHRAE ya se citan
en el texto. Antes de publicar:

- añadir NOM-008-ENER y NOM-020-ENER como normas, con organismo, edición, fecha,
  URL oficial, apéndice y fecha de consulta;
- completar DOI y número de fascículo cuando existan:
  - Xamán et al.: `10.1016/j.enbuild.2004.11.001`;
  - Hollands et al.: `10.1016/0017-9310(75)90179-9`;
  - Borbón et al.: `10.1612/inf.tecnol.4407it.10`;
  - Chow & Levermore: `10.1177/0143624407078642`;
  - Morillón et al.: `10.1016/j.solener.2003.11.008`;
- reemplazar `ASHRAE Handbook 1997` genérico por edición, capítulo, páginas y
  ecuación realmente usados, o por la fuente original de `RF`;
- verificar año/ISBN exactos de la sexta edición de Incropera usada;
- citar las fuentes de propiedades del aire;
- citar la especificación EPW y el proveedor concreto del archivo meteorológico;
- revisar que cada correlación se atribuya a la fuente donde aparece la ecuación,
  no solamente a una referencia que la cita.

Para citar el software, el artículo metodológico de 2016 y el ejecutable actual no
son el mismo objeto. La recomendación debe ser: citar el artículo para la genealogía
del método **y** una versión archivada del software para los resultados. Conviene
crear un DOI versionado —por ejemplo mediante Zenodo— y añadir a `CITATION.cff` la
fecha de lanzamiento, DOI y commit/release. `preferred-citation` no debería hacer
que la cita del software 0.2.1 desaparezca detrás del artículo del servicio de 2016.

## PyPI y mantenimiento de una única fuente

La página actual de [EnerHabitat 0.2.1 en PyPI](https://pypi.org/project/enerhabitat/)
se alimenta de `README.md`. PyPI no ofrece el mismo procesamiento matemático que
Quarto; las ecuaciones LaTeX pueden quedar como texto o renderizarse de forma
deficiente. Mantener en README un resumen sin depender de fórmulas complejas, o
usar expresiones de texto accesibles y enlazar la derivación completa.

Ahora hay duplicación sustancial entre README y las páginas Quarto. Para evitar que
errores de parámetros (como ocurrió con `h_i`) se corrijan en un sitio y persistan
en otro:

- establecer una página normativa única para símbolos, supuestos y parámetros;
- generar o incluir desde esa fuente los fragmentos repetidos;
- añadir una prueba que busque afirmaciones sensibles y tipos de retorno obsoletos;
- ejecutar ejemplos 1D y validar manifiestos 2D en CI;
- comprobar enlaces, referencias BibTeX y render de ecuaciones en cada release.

El tono promocional «we love [uv]» de README puede sustituirse por una instrucción
neutral. No es un problema técnico, pero mejora el registro de una publicación
científica.

## Arquitectura documental propuesta

1. **Alcance y límites:** qué resuelve, qué no resuelve y significado de las
   salidas.
2. **Entradas climáticas:** especificación EPW, día sintético, transposición solar,
   albedo, `RF`, fecha y zona horaria.
3. **Modelo físico 1D:** dominio, supuestos, signos, ecuaciones, setpoints y cargas.
4. **Modelo físico 2D:** dominio sólido, simetría, discretización geométrica,
   cavidades y rangos de correlaciones.
5. **Método numérico real:** ecuaciones discretas completas, acoplamiento,
   convergencia, estabilidad y convención temporal.
6. **Verificación y validación:** matriz por solver, datos, incertidumbre y
   resultados versionados.
7. **Guía de uso reproducible:** ejemplos con propiedades trazables y manifiestos.
8. **API:** referencia generada desde firmas/docstrings comprobadas, tipos y errores.
9. **Proveniencia y cita:** versiones, DOI, datos meteorológicos y licencias.

Conviene separar «Theory» de «Evidence». La validación no debe quedar como cuatro
viñetas al final de la página numérica; merece una página con casos, métricas,
figuras, datos descargables y límites de extrapolación.

## Plan de trabajo y criterios de aceptación

### Fase 0 — resolver antes del envío del artículo

- [ ] Revisar el promedio superficial 2D y la discretización de las fronteras.

**Criterio de salida:** ningún resultado se publica como convergido sin indicador y
residuales; cada ecuación de la documentación tiene correspondencia identificable
con el código; todas las diferencias de modelo frente a 2016 están incluidas en el
alcance de la evidencia o declaradas como no validadas.

### Fase 1 — exactitud documental y reproducibilidad

- [ ] Especificar por completo `meanDay()`, DNI/GHI/DHI, pvlib, albedo, `year` y
  zona horaria.
- [ ] Corregir tipos de retorno, firma de API, caché, `reset()` y diferencias 1D/2D.
- [ ] Publicar supuestos, signos, áreas de referencia y convenciones temporales.
- [ ] Añadir rangos de las correlaciones y fuentes de propiedades.
- [ ] Generar manifiestos y hashes de ejemplos; versionar el repositorio de
  validación y enlazar casos exactos.
- [ ] Completar bibliografía y `CITATION.cff`.

**Criterio de salida:** una persona externa puede reproducir cada número y saber
qué versión, datos, malla, tolerancia y criterio de éxito lo produjeron.

### Fase 2 — calidad editorial y mantenimiento

- [ ] Reducir duplicación README/Quarto.
- [ ] Verificar render de PyPI, enlaces, accesibilidad de figuras y glosario.
- [ ] Añadir una tabla de símbolos con unidades SI y una guía terminológica.
- [ ] Convertir benchmarks de rendimiento en experimentos reproducibles o retirar
  cifras.
- [ ] Añadir revisión automática de ejemplos y documentación en CI.

## Lista de comprobación final para una publicación científica

- [ ] El resumen no promete simulación del edificio ni consumo HVAC.
- [ ] Cada parámetro predeterminado tiene fuente o se declara elección del software.
- [ ] Cada correlación tiene geometría, longitud característica, propiedades y rango.
- [ ] Las ecuaciones continuas y discretas usan el mismo convenio de signos.
- [ ] Los límites de malla y la rasterización de geometría son cuantificables.
- [ ] La convergencia no usa promedios con signo y cubre todos los estados.
- [ ] Se informa cuando se alcanza un máximo sin converger.
- [ ] Hay cierre energético y estudios de independencia numérica.
- [ ] Regresión, verificación y validación se nombran correctamente.
- [ ] La evidencia de 2016 no se atribuye automáticamente a funciones nuevas.
- [ ] Los resultados 2D de techos tienen evidencia específica o se marcan como no
  validados experimentalmente.
- [ ] El día promedio, el año sintético y la zona horaria no pueden confundirse con
  una serie meteorológica real.
- [ ] Propiedades, EPW, resultados y código tienen procedencia, versión y hash.
- [ ] PyPI, README, sitio y artículo no contienen afirmaciones contradictorias.

## Fuentes externas principales consultadas

- [Barrios et al. (2016), artículo de referencia](https://doi.org/10.1016/j.solener.2015.12.017).
- [NOM-008-ENER-2001, PDF oficial](https://e.economia.gob.mx/wp-content/uploads/sites/29/PDF_Normas_Publicas/008ener.pdf).
- [Diccionario oficial del formato EPW de EnergyPlus](https://bigladdersoftware.com/epx/docs/22-1/auxiliary-programs/energyplus-weather-file-epw-data-dictionary.html).
- [pvlib `get_total_irradiance`, documentación oficial 0.13.1](https://pvlib-python.readthedocs.io/en/v0.13.1/reference/generated/pvlib.irradiance.get_total_irradiance.html).
- [Xamán et al. (2005)](https://doi.org/10.1016/j.enbuild.2004.11.001).
- [Hollands et al. (1975)](https://doi.org/10.1016/0017-9310(75)90179-9).
- [Borbón et al. (2010), texto completo](https://www.scielo.cl/pdf/infotec/v21n6/art17.pdf).
- [NIST, formulación de radiosidad para recintos grises difusos](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=841367).
- [EnerHabitat 0.2.1 en PyPI](https://pypi.org/project/enerhabitat/).
