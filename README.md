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

Paquetes usados: `amsmath`, `amssymb`, `amsthm`, `mathtools`, `enumitem`, `geometry`, `tikz`, `booktabs`, `hyperref`; `babel` (español) y `lmodern` se cargan solo si están instalados.

## Bibliografía principal

- I. R. Shafarevich, *Basic Algebraic Geometry*: libro de texto del curso.
- Palacios, texto sobre esquemas y haces coherentes (Capítulo 5: espectro y `Proj`).

La bibliografía de la monografía se redactó de memoria y conviene verificarla antes de citarla.

## Notas y puntos por verificar

- El campo `\author{}` está vacío en todos los documentos; complétese antes de difundirlos.
- El libro de Palacios afirma `φ*(V(𝔟)) = V(𝔟^c)`. Esto es falso en general: es cierto para la clausura de la imagen, o si `φ` es entera (contraejemplo: `k[x] → k[x]_x`). Está corregido en la guía y en las clases.
- Algunos resultados se usan sin demostrar y se indican como tales en el texto (Nullstellensatz, teorema de la base de Hilbert, teorema de Gauss sobre `A[x]`, `dim k[x₁,…,xₙ] = n`).

## Licencia

Por definir. Si el repositorio es público, se sugiere [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.es) para el texto.
