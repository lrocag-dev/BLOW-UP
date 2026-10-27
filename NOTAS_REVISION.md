# Notas de revisión

Observaciones detectadas al preparar el apéndice computacional (`apendice_computacional.tex`, `codigo/`).
**No se ha modificado el contenido matemático de ningún documento existente**: aquí solo se anotan los
puntos, con archivo, línea y justificación, para que el autor decida.

Convenciones: *Discrepancia* = algo que no coincide con lo comprobado; *Precisión* = el texto es correcto
pero admite una aclaración; *Documentación* = inconsistencia en la descripción del repositorio.

## Discrepancias en los documentos

No se encontraron errores matemáticos en lo que se verificó con cálculo (véase el apéndice, Secciones 2 y 4.1):
núcleo de `t -> (t², t³)`, potencias de `(x,y)`, función de Hilbert, cartas, resolución embebida de la cúspide
(`π*C = C̃ + 2E₁ + 3E₂ + 6E₃`, autointersecciones `-3,-2,-1`), fórmula de δ, paraguas de Whitney, tabla de
equilibrios, Proposición de equilibrios interiores, cartas `W₁` y `W₂` del infinito, los dos blow-ups y la
tabla de valores propios, pliegues para `Q = 2`, y las cantidades numéricas citadas (valores críticos
`0.3204, 0.3388, 1.3962`; raíces `0.266, 0.457, 0.660`; punto nilpotente `x*≈0.5751, y*≈0.3368, P≈0.6646, S≈0.3204`).

## Documentación del repositorio

1. **`README.md`, líneas 27 y 36 (y el orden de lectura)** — Se describe `guia_cap5_palacios` (9 págs.), pero ni `guia_cap5_palacios.tex`
   ni `.pdf` existen en el repositorio. Falta el archivo o sobra la entrada.
2. **`README.md`, línea 19** — `clase_blowup` figura con 23 págs.; el `clase_blowup.pdf` del repositorio tiene 14.
3. **`README.md`, línea 28** — `ejercicios_ideales` figura con 7 págs.; el `ejercicios_ideales.pdf` del
   repositorio tiene 5.
   (Posible causa de 2 y 3: PDFs compilados con otra versión de las fuentes. No se han tocado.)

## Precisiones (no son errores)

4. **`monografia_blowups.tex`, línea 1077, consecuencia (iii)** — Se dice que `S = 1` da una raíz y `S = 1.5`
   tres. Es cierto para las raíces *reales* de `h`, pero para `S = 1` (raíz `x ≈ 1.072`) y `S = 1.5`
   (raíces `≈ -1.064, -0.763, 1.183`) ninguna corresponde a un equilibrio del **primer cuadrante** (tienen
   `y = f(x) < 0` o `x < 0`). Conviene decirlo si el contexto es el primer cuadrante.
5. **`monografia_blowups.tex`, líneas 1099–1103 (tabla de los segundos blow-ups)** — Los valores propios
   dependen de por qué potencia se divide en cada carta. Por ejemplo, en la carta `ζ = vη` la dirección
   `v = PR·ζ` aparece en `η = 1/(PR)` con valores propios `{1/2, -1}`, no `{PR/2, -PR}` (que son los de la carta
   `v = ζξ`). Difieren por un reescalado del tiempo y el tipo (silla) es el mismo. Podría añadirse una frase.
6. **`monografia_blowups.tex`, líneas 919–926 (paraguas de Whitney)** — El texto es correcto: la transformada
   estricta es lisa y *no queda* singularidad. Si en otro lugar se espera «la singularidad restante», se
   comprobó con bases de Gröbner (`codigo/simbolico/07_paraguas_whitney.py`) que no existe; lo único especial
   es el punto pinzado (donde `X̃ → X` deja de ser inmersión).
7. **`monografia_blowups.tex`, línea 1077 (candidato Bogdanov–Takens)** — El valor `S ≈ 0.3204` coincide con
   el primer valor crítico de `Φ`; es decir, en ese punto `h` tiene una raíz *doble* (bifurcación silla-nodo).
   Es coherente (`Φ'(x*) = 0`), pero implica que la localización numérica del equilibrio es mal condicionada
   (error `~10⁻⁸` en `x*`): véase el apéndice, Sección 4.2.
