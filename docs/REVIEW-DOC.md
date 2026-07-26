# Revisión científica de la documentación de EnerHabitat

**Fecha de la revisión:** 24 de julio de 2026  
**Objeto revisado:** sitio Quarto en `docs/`, código al que remite la documentación,
pruebas, metadatos de publicación y artefactos precalculados del repositorio local  
**Base Git:** rama `main`, `99946d9`, más cambios locales sin confirmar  
**Estado del dictamen:** **revisión mayor antes de someter**  
**Revistas consideradas:** *SoftwareX* y *Journal of Building Performance
Simulation* (JBPS)

## 1. Dictamen ejecutivo

La documentación es extensa, clara y sustancialmente superior a la de la mayoría de
los paquetes científicos en una fase previa a publicación. La separación entre uso,
modelo 1D, modelo 2D, método numérico, supuestos y API es acertada. También son
fortalezas la declaración explícita del dominio de aplicación, la explicación de las
diferencias respecto de la herramienta de 2016, la exposición de los criterios de
convergencia y la inclusión de ejemplos reproducibles.

No obstante, **todavía no recomendaría presentar esta documentación como registro
científico definitivo del software**. Hay tres problemas que deben resolverse antes de
cualquier envío:

1. La sección “Validation” mezcla validación experimental o contra EnergyPlus de la
   herramienta de 2016 con verificación y pruebas de regresión del paquete actual.
   Debido a que la versión Python modificó entradas climáticas, malla, paso temporal,
   tolerancias y partes de la física de cavidades, la evidencia antigua no valida por
   sí sola la versión actual.
2. La instalación publicada no reproduce la documentación revisada: al 24 de julio
   de 2026, `pip install enerhabitat` instala 0.2.1 desde
   [PyPI](https://pypi.org/project/enerhabitat/), mientras el árbol local declara
   0.3.3 y otros metadatos declaran 0.3.0. Las diferencias entre esas versiones
   cambian resultados.
3. Las salidas 2D `Tso` y `Tsi`, documentadas como promedios superficiales, se
   calculan en los caminos de producción como `sum(T_surface)/(nx - 1)` aunque se
   suman `nx` valores y la malla usa `dx = X/nx`. El propio código califica esa
   convención como sesgada y usa `sum/nx` para decidir el coeficiente convectivo
   interior. Las gráficas y archivos precalculados que empleen esas columnas deben
   regenerarse después de corregirla.

Además, para JBPS falta una campaña documentada de independencia de malla y paso
temporal, sensibilidad numérica, incertidumbre y validación de la versión liberada.
Sin ella, el trabajo corre un riesgo alto de rechazo por insuficiencia de
verificación y validación. Para *SoftwareX*, la documentación se acerca más al tipo de
contribución esperado, pero la versión instalable, el archivo permanente del software
y el expediente reproducible de validación siguen siendo bloqueantes.

**Recomendación de destino:** en su estado actual, EnerHabitat tiene mejor encaje
potencial en *SoftwareX*, después de corregir los hallazgos P0 y P1. JBPS sería
razonable sólo si el artículo aporta una contribución metodológica de simulación más
allá de describir el paquete y presenta una campaña fuerte y actual de V&V,
incertidumbre y dominio de aplicabilidad.

## 2. Alcance y método de esta revisión

Se revisaron:

- la estructura, navegación, texto, ecuaciones, referencias, ejemplos y API del sitio;
- la correspondencia de las afirmaciones documentales con `src/enerhabitat/`;
- las pruebas disponibles y su relación con las afirmaciones de validación;
- los datos y resultados precalculados en `docs/data/results/`;
- `README.md`, `CHANGELOG.md`, `pyproject.toml` y `CITATION.cff`;
- la compilación local completa del sitio;
- las expectativas editoriales vigentes de ambas revistas.

La revisión se hizo sobre el árbol de trabajo local, que contiene cambios sin
confirmar; por ello este dictamen debe asociarse posteriormente a un *commit* y una
versión exactos. No se revisó un manuscrito científico actual: las observaciones sobre
las revistas se refieren a la documentación y al paquete que ésta pretende respaldar.

Las prioridades empleadas son:

- **P0 — bloqueante:** afecta la validez, los resultados o la reproducibilidad y debe
  resolverse antes del envío.
- **P1 — mayor:** probablemente generará una objeción sustantiva del revisor.
- **P2 — moderada/editorial:** no invalida el modelo, pero reduce claridad,
  mantenibilidad o calidad de presentación.

## 3. Fortalezas

### 3.1. Arquitectura de información

La navegación separa correctamente:

- preparación y configuración;
- ejemplos 1D y 2D;
- fundamentos físicos 1D y 2D;
- discretización, algoritmo y convergencia;
- supuestos y límites;
- referencia de API y citación.

La ruta de lectura funciona tanto para quien quiere ejecutar el programa como para
quien quiere auditar el modelo. La página `assumptions.qmd`, en particular, es una
buena decisión editorial: relaciona cada hipótesis con su consecuencia interpretativa
en lugar de ocultarla entre ecuaciones.

### 3.2. Transparencia del alcance físico

Se deja claro que se modela un solo componente opaco, sin ventanas, ventilación,
infiltración, ganancias internas, humedad, calor latente ni un sistema HVAC real.
También se explica que el modo acondicionado impone una temperatura y calcula carga
sensible ideal. Esta delimitación evita que el usuario confunda EnerHabitat con un
modelo integral de edificio.

Debe conservarse este nivel de franqueza. Sólo hace falta que la terminología de la
página inicial y los metadatos sea igual de precisa; véase P1-09.

### 3.3. Exposición del modelo y del método numérico

Son especialmente valiosos:

- la forma conservativa de las ecuaciones;
- la explicación del mapeo de capas subcelda en 1D;
- la distinción entre TDMA 1D y barridos iterativos 2D;
- los dos criterios no cancelables del lazo interior 2D;
- el cierre periódico día a día y los diagnósticos `day_error`,
  `inner_iterations`, `converged` y `energy_imbalance`;
- la explicación del nodo de aire interior y de los nodos de aire de cavidad;
- la distinción entre las correlaciones de Xamán y Hollands;
- la declaración de que Xamán se extrapola fuera del intervalo geométrico original;
- la tabla que compara la herramienta de 2016 con el paquete actual.

Estos elementos permiten una revisión científica real del programa y no sólo de su
interfaz.

### 3.4. Ejemplos y material de apoyo

Los ejemplos 1D son ejecutables y los 2D muestran geometría, campos y series
temporales. El inspector de secciones es útil para detectar errores de entrada antes
de un cálculo costoso. La compilación Quarto termina correctamente y genera las diez
páginas previstas sin errores de renderizado.

### 3.5. Evolución reciente del software

El árbol actual ya resuelve varias debilidades que habrían sido graves:

- cálculo de radiación por recinto en lugar de una suma par a par simplificada;
- criterios interiores no cancelables y diagnósticos de convergencia;
- mapeo 1D sensible a interfaces y conservación de masa térmica;
- selección de `hi` para techos según dirección del flujo;
- página explícita de supuestos y diferencias respecto de 2016;
- pruebas de reducción 2D→1D, balance de energía y regresión.

La revisión que sigue no niega ese avance. Precisamente porque el paquete cambió de
manera material, la evidencia científica y la versión publicada deben actualizarse a
la misma velocidad.

## 4. Adecuación a las revistas

| Criterio | *SoftwareX* | JBPS | Estado observado |
|---|---|---|---|
| Software accesible, inspeccionable y reutilizable | Central | Importante | Código abierto y licencia MIT: favorable |
| Versión exacta y archivable | Central | Central para reproducibilidad | No resuelto: PyPI 0.2.1, árbol 0.3.3, CFF/tag 0.3.0 |
| Material de apoyo y caso ejecutable | Central | Importante | Buena base, pero los resultados 2D carecen de manifiesto |
| Verificación y validación del código liberado | Importante | Central | Se confunden evidencias de 2016, regresión y validación actual |
| Independencia numérica y cuantificación de incertidumbre | Deseable | Muy importante | Insuficiente, especialmente en 2D |
| Novedad metodológica en simulación de edificios | Secundaria | Central | No puede inferirse sólo de la documentación |
| Reproducibilidad de código, datos y entorno | Central | Explícitamente promovida | Aún incompleta |

La página oficial de [*SoftwareX*](https://www.sciencedirect.com/journal/softwarex)
señala que el software debe estar públicamente disponible para inspección, validación
y reutilización, acompañado por su distribución de código y material de soporte.
La página oficial de
[JBPS](https://www.tandfonline.com/journals/tbps20) incluye expresamente calidad,
usabilidad, validación, verificación, cuantificación de incertidumbre y flujos de
software entre sus temas. Su editorial sobre
[reproducibilidad en simulación de edificios](https://www.tandfonline.com/doi/abs/10.1080/19401493.2024.2441385)
subraya la necesidad de código, datos, entorno y métodos detallados.

**Conclusión editorial:** la documentación puede convertirse en un buen paquete de
soporte para *SoftwareX*. Para JBPS no basta con ampliar el texto; hace falta producir
evidencia científica adicional.

## 5. Hallazgos P0 — bloqueantes

### P0-01. La “validación” no corresponde inequívocamente al software actual

`numerics.qmd` afirma que el registro de validación combina la validación publicada
con pruebas de regresión del paquete. Después presenta, en una sola lista:

- comparación 1D contra EnergyPlus publicada en 2016;
- comparación de muro hueco contra caja caliente publicada;
- reducción 2D→1D del paquete;
- reproducción de salidas del código C mediante *golden masters*.

Estas evidencias responden a preguntas distintas:

- **validación:** ¿el modelo representa adecuadamente el fenómeno físico?;
- **verificación de solución:** ¿la discretización converge a la solución matemática?;
- **verificación de código:** ¿la implementación resuelve las ecuaciones previstas?;
- **regresión:** ¿el código conserva una salida histórica?;
- **equivalencia:** ¿dos implementaciones producen resultados concordantes?

Un *golden master* puede detectar cambios, pero no demuestra que la salida sea
físicamente correcta. La reducción 2D→1D es una prueba de consistencia, no una
validación externa. Y la validación de la herramienta de 2016 no se transfiere
automáticamente a la versión actual porque ésta emplea, entre otras diferencias:

- EPW y transposición solar con pvlib en lugar de la entrada climática histórica;
- `dt = 10 s` en lugar de 1 s;
- otras mallas y tolerancias;
- nuevo mapeo de capas;
- tratamiento actualizado de radiación;
- otra reducción de la correlación de Xamán;
- selección de `hi` por flujo y ajustes del factor radiativo de onda larga.

`CHANGELOG.md` menciona una campaña independiente en otro repositorio, pero la
documentación no proporciona un enlace público, identificador de versión, datos,
scripts ni resultados auditables.

**Acción exigida:**

1. Crear una página “Verification and validation” independiente de “Numerical
   method”.
2. Separar explícitamente las cinco categorías anteriores.
3. Publicar y enlazar un repositorio o archivo versionado de validación.
4. Ejecutar la campaña con la misma versión que se cite y distribuya.
5. Incluir por caso: entrada bruta, preprocesamiento, versión, configuración completa,
   referencia, métrica, incertidumbre, criterio de aceptación y resultado.
6. Mientras no exista esa evidencia, sustituir “the package is validated” por una
   formulación precisa: la herramienta de 2016 fue validada; la versión actual tiene
   determinadas pruebas de verificación; su validación independiente está pendiente o
   se encuentra en el recurso enlazado.

Una matriz mínima debería verse así:

| Caso | Propósito | Modelo/modo | Referencia independiente | Métricas mínimas |
|---|---|---|---|---|
| Paredes multicapa | Validación | 1D libre y AC | EnergyPlus o benchmark analítico | amplitud, fase, `Ti`, calor diario |
| Muro relleno homogéneo | Verificación cruzada | 2D→1D | solver 1D refinado | normas L∞/L2, energía |
| Bloque hueco | Validación | 2D aire | caja caliente, datos brutos | flujo, U/R, incertidumbre |
| Bovedilla | Validación | 2D aire/solid | experimento o benchmark independiente | flujo y temperaturas |
| Todos los anteriores | Verificación numérica | todos | solución refinada | malla, `dt`, tolerancias |

### P0-02. La instalación indicada no reproduce la documentación

`index.qmd` y `README.md` indican:

```bash
pip install enerhabitat
```

Al momento de esta revisión, ese comando instala 0.2.1 desde PyPI. El árbol local
declara 0.3.3 en `pyproject.toml`; `CITATION.cff`, el valor de respaldo de
`src/enerhabitat/__init__.py` y la etiqueta disponible declaran 0.3.0; el
*changelog* visible en `about.qmd` termina en 0.2.1. No existe por tanto un objeto
único que conecte documentación, código, resultados y cita.

Este problema es científico, no sólo editorial: `CHANGELOG.md` reconoce cambios que
alteran resultados entre 0.2.1 y 0.3.x.

**Acción exigida:**

- decidir la versión que respaldará el artículo;
- sincronizar `pyproject.toml`, `CITATION.cff`, `__version__`, etiqueta Git, PyPI,
  documentación y changelog;
- publicar esa versión antes de pedir a los revisores que ejecuten los ejemplos;
- fijar el comando del artículo y del sitio, por ejemplo
  `pip install enerhabitat==X.Y.Z`;
- mostrar versión y *commit* de construcción en el pie o encabezado del sitio;
- crear un archivo permanente de la liberación, preferentemente con DOI;
- citar tanto el artículo metodológico de 2016 como el objeto de software exacto.

### P0-03. `Tso` y `Tsi` 2D no son los promedios que describe la API

`api.qmd` y `usage-2d.qmd` llaman a `Tso` y `Tsi` “outer/inner-surface mean”.
Sin embargo, los caminos de producción en `src/enerhabitat/ehtools2d.py` suman los
`nx` nodos de la superficie y dividen entre `nx - 1`. Esto ocurre en las variantes
libre/AC, sólida/con aire y muro/losa. La malla se construye con `dx = X/nx`.

El problema no es hipotético: comentarios del mismo archivo describen `sum/nx` como
el promedio no sesgado para decidir `hi` y califican la convención reportada
`/(nx - 1)` como sesgada. Con `nx = 80`, la serie publicada se multiplica por
`80/79`. Aplicar ese factor directamente a temperaturas Celsius produce un error de
alrededor de 0.3 °C para valores próximos a 25 °C y, además, depende del origen de la
escala. Si se pretendiera una integración trapezoidal sobre nodos de frontera,
deberían existir pesos de media celda en los extremos; no existen en la suma actual.

**Acción exigida:**

1. Definir formalmente dónde se ubican los grados de libertad y qué integral
   representa la temperatura superficial.
2. Corregir el promedio de producción; con la convención de volúmenes actual,
   `sum/nx` es la opción consistente.
3. Añadir una prueba invariante: un campo superficial uniforme de valor `C` debe
   producir exactamente `C` para cualquier `nx` y en Kelvin o Celsius.
4. Regenerar todos los CSV, NPY, figuras congeladas y números del manuscrito que usen
   `Tso` o `Tsi`.
5. Mantener la convención histórica sólo dentro de pruebas de fidelidad al código C,
   claramente separada de las salidas científicas.

### P0-04. Debe cerrarse la formulación de los volúmenes de frontera

En 1D se documentan celdas centradas de ancho completo `dx = L/Nx`, con capacidad
`rho*c*dx`. En el ensamblaje se añade directamente `h_o` o `h_i` al volumen de
frontera. En 2D se sigue una operación análoga con `h*dx`. No se documenta si el
grado de libertad exterior está en el centro de la celda o en la superficie física.

Si está en el centro, el intercambio entre ese centro y el ambiente incluye la
resistencia de media celda y la película:

```text
R_boundary = (dx/2)/k + 1/h
```

Si está en la superficie, debe justificarse la capacidad y las distancias a las
caras vecinas de ese volumen de frontera. La documentación afirma que los
coeficientes entran “exactly as the boundary conditions prescribe”, pero no aporta
el diagrama ni la derivación necesarios para verificarlo.

No afirmo aquí que la solución correcta sea necesariamente una sola modificación de
coeficiente; afirmo que la convención actual es ambigua y que sus ecuaciones
discretas no permiten auditar la consistencia geométrica.

**Acción exigida:**

- añadir un esquema de centros, caras, anchos y medios volúmenes;
- derivar los coeficientes de los nodos exterior e interior;
- verificar la implementación con una pared homogénea estacionaria cuya resistencia
  exacta sea `1/ho + L/k + 1/hi`;
- realizar refinamiento de malla para comprobar que el flujo y las temperaturas
  superficiales convergen a ese valor;
- corregir código y resultados si el ensayo revela la omisión de media celda.

## 6. Hallazgos P1 — mayores

### P1-01. Falta independencia numérica sistemática

Existe una prueba 1D que compara `Nx = 200` y `Nx = 800` y exige una diferencia máxima
de `Ti` inferior a 0.05 °C, pero ese resultado no se presenta en la documentación.
No se encontró un estudio equivalente y documentado para:

- malla 2D en ambas direcciones;
- `dt`;
- `tol_inner`, `tol_day`, `max_inner` y `max_days`;
- rasterización de paredes, nervaduras y cavidades;
- modos libre y AC;
- geometrías `HollowBlock` y `Slab`, con relleno de aire y sólido.

La estabilidad del esquema implícito de conducción no demuestra exactitud. Para JBPS
este punto es bloqueante.

**Acción recomendada:** presentar tablas de refinamiento para casos representativos,
con costo computacional, `Ti` extrema y amplitud, desfase, energía diaria,
temperaturas de superficie/cavidad, normas de error y razón para elegir los valores
por omisión.

### P1-02. El dominio de estabilidad de los nodos de cavidad no se controla

La documentación reconoce que los nodos de aire de cavidad se actualizan
explícitamente y que `lambda_h*dt` puede acercarse a 1 en el ejemplo de losa. La API,
sin embargo, acepta geometrías que pueden producir áreas pequeñas, perímetros grandes
o coeficientes convectivos altos, sin calcular ni imponer el límite de estabilidad.

**Acción recomendada:**

- calcular y registrar el máximo `lambda_h*dt` durante la simulación;
- advertir o rechazar valores fuera del dominio establecido;
- aplicar subpasos o una actualización implícita si se quiere admitir geometrías más
  generales;
- incluir el rango de Rayleigh y la razón de aspecto observados, no sólo las fórmulas
  de correlación.

### P1-03. La geometría y la emisividad no tienen un contrato de entrada suficiente

La API enumera claves geométricas, pero no especifica de forma completa:

- positividad y relaciones entre espesores;
- límites de `n_cavities`;
- condición `topping_cap <= topping`;
- resolución mínima de paredes y cavidades en la malla;
- intervalo físico `0 < emissivity <= 1`;
- comportamiento ante dimensiones que desaparecen al rasterizarse;
- rango de validez de las correlaciones.

La validación previa al cálculo comprueba capas, orientación y configuración general,
pero no todas esas invariantes.

**Acción recomendada:** añadir una tabla de dominio admisible, validar cada relación
y reportar dimensiones efectivas discretizadas junto con el error respecto de la
geometría solicitada.

### P1-04. Se afirma soporte para emisividades por superficie que la API no implementa

`model-2d.qmd` afirma que “per-surface emissivities are supported”. En producción,
`_transfer_factors(..., emissivity)` convierte el argumento a un único `float`; los
constructores y la API también documentan un escalar uniforme.

**Acción recomendada:** eliminar la frase o implementar y probar un vector de
emisividades con una convención de orden inequívoca. La documentación no debe atribuir
a la formulación matemática una capacidad ausente en la interfaz y el código.

### P1-05. Algunas afirmaciones de convergencia son más fuertes que la métrica

La página numérica dice que “every persisted state stops changing”, pero la parte del
campo se mide con el promedio nodal de cambios absolutos, no con la norma máxima.
Una región pequeña puede cambiar más que la tolerancia y aun así satisfacer el
promedio. Los nodos de aire sí entran mediante máximos escalares.

Asimismo, la estimación `8640 * 1e-8 ≈ 1e-4 °C` no constituye, sin una demostración de
estabilidad o contracción, una cota rigurosa del error acumulado del día.

**Acción recomendada:** describir literalmente la norma empleada; reportar también
el máximo nodal; justificar empíricamente la tolerancia; y definir un umbral recomendado
para `energy_imbalance`. Si se desea sostener “cada estado”, usar una norma máxima.

### P1-06. La convención temporal mezcla estados de inicio y fin de paso

La API reconoce que `Tso` se guarda antes de resolver el paso y `Ti`/`Tsi` después,
pero todas las variables reciben el mismo índice de tiempo. El 1D también aplica el
forzamiento del instante etiquetado, avanza el estado y guarda la respuesta bajo ese
mismo índice.

Esto es relevante para desfases, picos, comparación experimental y validación contra
otro programa.

**Acción recomendada:** adoptar una sola convención —estado al inicio, al final o en
el centro del intervalo—, etiquetar los resultados en el instante físico
correspondiente y documentar la cuadratura de energía. Si se conserva la convención
histórica, proporcionar columnas o coordenadas temporales distintas y cuantificar el
desfase de `dt`.

### P1-07. El día promedio necesita una especificación reproducible de los intervalos EPW

La explicación climática ha mejorado: identifica GHI/DNI/DHI, zona horaria y que EPW
acumula radiación durante el intervalo anterior. Falta cerrar dos detalles:

- se interpolan las muestras horarias de 00:00 a 23:00 y después se prolonga el valor
  de las 23:00 hasta el final del día, sin cierre cíclico hacia 00:00;
- los valores de energía horaria de EPW se tratan como muestras de potencia situadas
  en la hora, por lo que debe demostrarse si la interpolación conserva el total diario
  pretendido.

**Acción recomendada:** documentar el algoritmo con pseudocódigo, añadir una figura de
las últimas horas del día y una prueba de conservación de la irradiación diaria.
Indicar explícitamente si se desea una señal periódica continua en el cambio de día.

### P1-08. Los resultados 2D precalculados no son artefactos autodescriptivos

`docs/data/results/summary.json` guarda tiempo, número de días, dimensiones y energías.
No guarda:

- versión y *commit*;
- versión de Python y dependencias;
- sistema operativo, procesador y contexto de JIT;
- `nx`, `ny`, `dt` y tolerancias;
- `converged`, `day_error`, `inner_iterations` y `energy_imbalance`;
- huellas de EPW, materiales y configuración;
- fecha y comando de generación.

Además, `freeze: auto` no garantiza por sí mismo que cambios en CSV/NPY o datos
auxiliares invaliden correctamente todas las páginas.

**Acción recomendada:** generar un manifiesto por corrida y hacer que el render falle
si no coincide con el software y los archivos actuales. Por ejemplo:

```yaml
case: slab_free
enerhabitat_version: X.Y.Z
git_commit: "<sha>"
python: "3.x.y"
dependencies_lock_sha256: "<sha256>"
inputs:
  epw_sha256: "<sha256>"
  materials_sha256: "<sha256>"
solver:
  nx: 80
  ny: 160
  dt_s: 10
  tol_inner: 1.0e-8
  tol_day: 5.0e-4
diagnostics:
  converged: true
  days: 6
  day_error_C: ...
  max_inner_iterations: ...
  energy_imbalance: ...
```

### P1-09. “Demanda de aire acondicionado” sobrestima el alcance de la salida

La portada dice que el paquete produce “air-conditioning energy demands”. En realidad
calcula carga térmica sensible ideal, por unidad de área, a través de un único
componente opaco, con temperatura interior prescrita. No incluye rendimiento de
equipo, energía eléctrica, distribución, ventilación, cargas internas, humedad ni
interacciones del edificio completo.

**Acción recomendada:** emplear desde el título y la descripción breve “ideal sensible
heating/cooling load through one opaque component” y reservar “building demand” o
“consumption” para modelos que incluyan el resto del sistema.

### P1-10. Se promete una matriz de ejemplos que no está completa

`README.md` afirma que las páginas contienen la matriz completa
1D/2D × libre/AC y los ejemplos de muro y losa. `run_examples.py` y los resultados
incluyen:

- `hollow_free`;
- `hollow_ac`;
- `slab_free`.

No existe `slab_ac`.

**Acción recomendada:** añadir el caso o restringir la afirmación. La documentación
debe contener al menos un ejemplo ejecutado de cada modo cuya equivalencia API se
promete.

### P1-11. Las afirmaciones de rendimiento no son reproducibles

La documentación menciona 10–20 minutos para un cálculo 2D, pero el propio resumen
registra aproximadamente 25.8 minutos para `hollow_free`. También declara aceleraciones
de ≈6× con 8 procesos y ≈10× con 16 sin protocolo, dispersión, hardware detallado ni
datos. `run_examples.py` ejecuta tres casos simultáneos, por lo que sus tiempos sufren
contención e incluyen compilación JIT.

La frase que dice que aumentar `dt` aporta “almost no speed” porque aumentan las
iteraciones interiores tampoco es general: el solver 1D no tiene iteraciones
interiores.

**Acción recomendada:** eliminar cifras no respaldadas o publicar un benchmark con
versiones, hardware, calentamiento JIT, número de repeticiones, mediana/variación y
casos separados. Diferenciar tiempo de primera ejecución y tiempo caliente.

### P1-12. Procedencia insuficiente de clima y materiales

El EPW incluido identifica en su cabecera una fuente TMYx, pero la documentación no
registra de forma completa URL de origen, fecha de descarga, licencia o condiciones de
redistribución y huella criptográfica. Los materiales se presentan correctamente como
ilustrativos, pero no incluyen fuentes, temperatura de referencia ni incertidumbre.

**Acción recomendada:** distinguir datos tutoriales de datos de validación; añadir un
archivo de procedencia y licencias; registrar SHA-256; y no usar propiedades
“típicas” como evidencia del artículo sin fuente y análisis de sensibilidad.

### P1-13. La cadena de pruebas no es directamente reproducible por un revisor

El proyecto no declara un grupo de dependencias de pruebas; `uv run pytest -q` falla
porque `pytest` no está instalado. Algunas pruebas de fidelidad se omiten cuando no se
dispone de las fuentes C históricas, y un *golden* de cavidad se marca como no
generado. Tampoco se encontró un flujo de integración continua en el árbol revisado.

**Acción recomendada:**

- añadir un grupo `test` y un comando único documentado;
- configurar CI para las versiones de Python compatibles;
- distinguir en el reporte las pruebas siempre ejecutables de las históricas;
- archivar legalmente los insumos C necesarios o reemplazarlos con vectores de prueba
  autosuficientes;
- publicar el reporte de pruebas asociado a la liberación.

## 7. Hallazgos P2 — moderados y editoriales

### P2-01. Citación y bibliografía

- La referencia de Borbón et al. sólo incluye una URL al PDF; conviene añadir DOI si
  se confirma en la fuente editorial.
- NOM-020 carece de URL y localización exacta del apartado empleado.
- La función por tramos que genera `DeltaTn` necesita una cita primaria con página,
  tabla o ecuación; citar un atlas de forma general no basta para reconstruir los
  umbrales.
- En `model-1d.qmd` se enuncia continuidad de temperatura y flujo, pero se muestra de
  forma inequívoca sólo la segunda. Añádase la igualdad de temperaturas a ambos lados
  de la interfaz.
- La página “About” debe citar por separado el método de 2016 y la liberación actual
  del software.

### P2-02. Accesibilidad y referencias a figuras

El HTML generado contiene figuras sin atributo alternativo. Varias salidas de código
carecen de leyenda numerada y referencia cruzada desde el texto.

**Acción recomendada:** añadir `fig-cap`, `fig-alt` e identificadores Quarto a todas
las figuras; describir la conclusión visual en el texto; comprobar contraste y
navegación por teclado.

### P2-03. Riesgo de deriva de la API manual

La referencia de API se mantiene a mano y duplica firmas, valores por omisión y
descripciones del código. La contradicción sobre emisividades muestra el riesgo.

**Acción recomendada:** generar al menos las firmas desde los objetos o añadir pruebas
que comparen documentación y API (`inspect.signature`, valores por omisión y nombres
de atributos).

### P2-04. Información de versión y vigencia del sitio

El sitio no muestra versión documentada, *commit*, fecha de actualización ni selector
de versiones. El changelog de “About” está atrasado respecto del repositorio.

**Acción recomendada:** añadir un encabezado o pie global y conservar documentación de
las versiones publicadas, evitando que una URL estable describa silenciosamente una
versión no liberada.

### P2-05. Guía de contribución y gobierno del proyecto

Para una publicación de software conviene incluir:

- instalación de desarrollo;
- estilo, pruebas y construcción de documentación;
- procedimiento para reportar errores científicos;
- política de versiones y compatibilidad;
- mantenimiento previsto y responsables;
- política de archivo de datos y liberaciones.

### P2-06. Terminología y consistencia

Revisar de manera global:

- “average day” frente a “representative day” y su definición estadística;
- “energy demand”, “load”, “heat transfer” y “consumption”;
- energía diaria en J/m² frente a unidades más legibles como kWh/m²·día;
- “node” frente a “control volume” y “surface node”;
- “validated”, “verified”, “regression-tested” y “legacy-compatible”.

La consistencia terminológica es especialmente importante porque cada término tiene
un significado técnico distinto en V&V.

## 8. Plan priorizado antes del envío

### Fase A — congelar el objeto científico

1. Elegir una versión única.
2. Corregir `Tso`/`Tsi` y resolver la formulación de frontera.
3. Ejecutar nuevamente todas las pruebas y regenerar resultados.
4. Sincronizar PyPI, etiqueta, CFF, changelog y documentación.
5. Crear liberación archivada y registrar su DOI o identificador permanente.

### Fase B — producir evidencia

1. Ejecutar independencia de malla, `dt` y tolerancias.
2. Publicar la validación actual contra EnergyPlus y caja caliente.
3. Añadir un caso independiente para `Slab` o delimitar expresamente que no está
   validado.
4. Documentar incertidumbre de datos y sensibilidad de propiedades/correlaciones.
5. Publicar matrices, datos brutos, scripts y criterios de aceptación.

### Fase C — reproducibilidad

1. Añadir entorno de prueba y CI.
2. Crear manifiestos de corridas con huellas criptográficas.
3. Fijar versiones en los comandos del artículo.
4. Ejecutar el flujo completo desde un clon limpio.
5. Depositar código, datos y resultados de la versión presentada.

### Fase D — revisión editorial

1. Separar la página de V&V del método numérico.
2. Ajustar alcance y terminología de cargas.
3. Corregir promesas de cobertura y rendimiento.
4. Completar citas, texto alternativo y referencias cruzadas.
5. Hacer una revisión de idioma científico consistente con el manuscrito.

## 9. Criterios mínimos de aceptación de una revisión

Yo consideraría resuelta esta revisión cuando se pueda responder “sí” a todo lo
siguiente:

- [ ] El comando de instalación obtiene exactamente la versión documentada.
- [ ] Código, etiqueta, PyPI, CFF, changelog, sitio y DOI muestran la misma versión.
- [ ] `Tso` y `Tsi` pasan una prueba de campo uniforme y se regeneraron los artefactos.
- [ ] La discretización de frontera reproduce la resistencia estacionaria analítica.
- [ ] Existe una matriz pública de verificación y validación de la versión liberada.
- [ ] Validación, verificación, regresión y fidelidad histórica están separadas.
- [ ] Se publican estudios de malla, tiempo y tolerancias para 1D y 2D.
- [ ] El solver controla o reporta el parámetro de estabilidad de cada cavidad.
- [ ] Toda geometría y emisividad se valida contra un dominio documentado.
- [ ] Cada resultado precalculado tiene versión, configuración, huellas y diagnósticos.
- [ ] Un clon limpio puede construir el sitio y ejecutar las pruebas con un comando.
- [ ] Las afirmaciones de rendimiento proceden de un benchmark reproducible.
- [ ] El alcance se describe como carga sensible ideal de un componente opaco.
- [ ] Las figuras son accesibles, numeradas y referidas en el texto.
- [ ] La documentación cita el artículo metodológico y la versión exacta del software.

## 10. Recomendación final por revista

### *SoftwareX*

**Revisión mayor, con perspectiva favorable después de corregirla.** La documentación
ya explica el propósito, arquitectura de uso, física y ejemplos con una profundidad
adecuada para acompañar un artículo corto de software. Antes de someter son
imprescindibles una liberación instalable y archivable, la corrección de salidas 2D,
un flujo reproducible y un expediente de validación atribuible a esa liberación.

### *Journal of Building Performance Simulation*

**No someter todavía sobre la base de esta documentación.** El sitio es un buen
material suplementario, pero no sustituye una contribución científica con V&V y
cuantificación de incertidumbre. Para hacer viable JBPS se requiere, como mínimo, una
campaña actual y convincente de validación, convergencia espacial/temporal, sensibilidad
de los parámetros de cavidad, incertidumbre y delimitación del dominio de validez.
También debe quedar explícita la novedad respecto de Barrios et al. (2016).

En síntesis: **la documentación es prometedora y técnicamente seria, pero la cadena
“versión exacta → ecuaciones → implementación → pruebas → validación → resultados
publicados” aún no está cerrada**. Cerrar esa cadena debe preceder al pulido estilístico
y a la selección definitiva de revista.
