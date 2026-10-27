# Blow-ups, blow-downs y resolución de singularidades

Proyecto de investigación en geometría algebraica a nivel doctoral: material de estudio, monografía y colecciones de ejercicios resueltos sobre **blow-ups, blow-downs y resolución de singularidades**, desarrollados con anillos graduados y morfismos de esquemas.

Todo el contenido está en español y escrito en LaTeX.

## Contenido

### Monografía

| Archivo | Páginas | Descripción |
|---|---|---|
| `monografia_blowups` | 52 | Monografía principal. Geometría algebraica pura, con las demostraciones centrales completas. Cubre la propiedad universal del blow-up, el teorema del centro liso (fórmula `K = π*K + (c−1)E`), la fórmula de δ, resolución embebida de curvas, la resolución de las singularidades `A_n`, el paraguas de Whitney y la aplicación a un sistema dinámico polinomial plano. Cita sin demostrar los teoremas de Castelnuovo, Zariski y Hironaka. |

### Clases

| Archivo | Páginas | Descripción |
|---|---|---|
| `clase_blowup` | 23 | Clase autocontenida de dos horas sobre blow-ups. Es el origen de la monografía. |
| `clase_ideales` | 10 | Sesión de 2 horas sobre ideales: ideal propio, generado, primo, maximal, radical, correspondencia y contracción. |
| `clase_cap5_espectro` | 12 | Clase sobre el Capítulo 5 de Palacios: `Spec`, topología de Zariski, `φ*`, localización y `Proj`. |

### Apéndice computacional

| Archivo | Páginas | Descripción |
|---|---|---|
| `apendice_computacional` | 13 | Experimentos computacionales que verifican e ilustran los ejemplos centrales: bases de Gröbner en aritmética exacta (núcleo de `t ↦ (t², t³)`, potencias de `(x,y)`, función de Hilbert, cartas del blow-up, resolución embebida de la cúspide, δ, paraguas de Whitney), figuras y una simulación numérica de las cartas del blow-up del sistema dinámico de la monografía. Tablas de *esperado* frente a *calculado* y sección «Alcance y límites». **Los cálculos no sustituyen demostraciones.** Es autocontenido: no depende de que existan los otros `.tex` (las referencias a la monografía son textuales). |

### Código

| Carpeta | Descripción |
|---|---|
| `codigo/` | Código reproducible del apéndice: `README.md` (cómo ejecutar, versiones), `requirements.txt` (versiones fijadas), `ejecutar_todo.sh`, `registro_ejecucion.txt` (salida de referencia). |
| `codigo/simbolico/` | Cálculo simbólico exacto (sympy; contraste opcional con Singular). |
| `codigo/figuras/` | Scripts de las figuras; salida (PDF vectorial) en `codigo/figuras/figuras_salida/`. |
| `codigo/dinamica/` | Verificación simbólica del sistema dinámico, equilibrios y simulación numérica en las cartas del blow-up (evidencia numérica). |

### Guía y ejercicios

| Archivo | Páginas | Descripción |
|---|---|---|
| `guia_cap5_palacios` | 9 | Guía del Capítulo 5 de Palacios: ejemplos del libro explicados, ejemplos nuevos, tabla de espectros y lista de erratas. |
| `ejercicios_ideales` | 7 | 25 ejercicios resueltos: ideales, primos, maximales, principales, generados, radicales, cocientes, cuerpo de fracciones, DIP y DFU. |
| `ejercicios_espectro` | 7 | 20 ejercicios resueltos: espectro de un anillo, topología de Zariski, `D(f)`, `φ*`, localización, irreducibilidad, `Proj` y cartas del blow-up. |

Cada documento existe como `.tex` (fuente) y `.pdf` (compilado).

## Orden de lectura sugerido

1. `clase_ideales` y `ejercicios_ideales`: vocabulario algebraico básico.
2. `clase_cap5_espectro`, `guia_cap5_palacios` y `ejercicios_espectro`: el espectro, la topología de Zariski y `Proj`.
3. `clase_blowup`: primera visión del blow-up.
4. `monografia_blowups`: desarrollo completo.
5. `apendice_computacional` (y `codigo/`): verificación con cálculo de los ejemplos centrales.

## Compilación

Se necesita una distribución de TeX con `pdflatex` (TeX Live o MiKTeX). Cada documento es autocontenido y se compila de forma independiente; conviene ejecutar dos pasadas para resolver las referencias:

```bash
pdflatex monografia_blowups.tex
pdflatex monografia_blowups.tex
```

O todos a la vez:

```bash
for f in *.tex; do pdflatex -interaction=nonstopmode "$f" && pdflatex -interaction=nonstopmode "$f"; done
```

Paquetes usados: `amsmath`, `amssymb`, `amsthm`, `mathtools`, `enumitem`, `geometry`, `tikz`, `booktabs`, `hyperref`; `babel` (español) y `lmodern` se cargan solo si están instalados. El apéndice computacional usa además `graphicx` y `fancyvrb` (y `tikz-cd`, como las clases).

## Cómo reproducir los resultados del apéndice computacional

Se necesita Python (probado con 3.13) y, para compilar el apéndice, `pdflatex`. Desde la raíz del repositorio:

```bash
# 1) Todos los cálculos (crea codigo/.venv, instala las versiones fijadas y ejecuta todos los scripts)
bash codigo/ejecutar_todo.sh

# 2) El documento (dos pasadas; usa las figuras y la tabla que producen los scripts)
pdflatex apendice_computacional.tex
pdflatex apendice_computacional.tex
```

Cada script termina con código 0 solo si pasan todas sus comprobaciones. Para ejecutar uno suelto:

```bash
python3 -m venv codigo/.venv
codigo/.venv/bin/pip install -r codigo/requirements.txt
codigo/.venv/bin/python codigo/simbolico/05_resolucion_cuspide.py
```

Detalles (versiones, semillas, tolerancias, qué produce cada script): [`codigo/README.md`](codigo/README.md).
Observaciones sobre los documentos existentes detectadas durante este trabajo: [`NOTAS_REVISION.md`](NOTAS_REVISION.md).

## Bibliografía principal

- I. R. Shafarevich, *Basic Algebraic Geometry*: libro de texto del curso.
- Palacios, texto sobre esquemas y haces coherentes (Capítulo 5: espectro y `Proj`).

La bibliografía de la monografía se redactó de memoria y conviene verificarla antes de citarla.

## Notas y puntos por verificar

- El campo `\author{}` está vacío en todos los documentos; complétese antes de difundirlos.
- El libro de Palacios afirma `φ*(V(𝔟)) = V(𝔟^c)`. Esto es falso en general: es cierto para la clausura de la imagen, o si `φ` es entera (contraejemplo: `k[x] → k[x]_x`). Está corregido en la guía y en las clases.
- Algunos resultados se usan sin demostrar y se indican como tales en el texto (Nullstellensatz, teorema de la base de Hilbert, teorema de Gauss sobre `A[x]`, `dim k[x₁,…,xₙ] = n`).
- Los cálculos de `codigo/` son verificación e ilustración, no demostraciones; los numéricos (punto flotante) no deciden igualdades. Las inconsistencias detectadas en los documentos existentes (sin corregirlas) están en `NOTAS_REVISION.md`.

## Licencia

Por definir. Si el repositorio es público, se sugiere [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.es) para el texto.
