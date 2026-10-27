# Código del apéndice computacional

Código reproducible detrás de [`../apendice_computacional.tex`](../apendice_computacional.tex). Verifica e ilustra
con cálculo los ejemplos centrales de la monografía. **Los cálculos no sustituyen demostraciones.**

- Todo lo algebraico se hace en **aritmética exacta** (sympy; racionales y cuerpos finitos $\mathbb F_p$). Nunca flotantes.
- Lo numérico (`dinamica/`) es **punto flotante** y se declara como evidencia numérica.
- Todas las rutas son relativas al propio script; no hay descargas fuera de PyPI (`pip`).

## Ejecución rápida

Desde la raíz del repositorio:

```bash
bash codigo/ejecutar_todo.sh
```

Crea `codigo/.venv`, instala `requirements.txt`, ejecuta todos los scripts y escribe el registro en
`codigo/registro_ejecucion.txt`. Falla (código de salida ≠ 0) si no pasa alguna comprobación.

### Ejecución manual

```bash
python3 -m venv codigo/.venv
codigo/.venv/bin/pip install -r codigo/requirements.txt
codigo/.venv/bin/python codigo/simbolico/01_nucleo_cuspide.py      # y así con cada script
```

## Versiones

| | |
|---|---|
| Python | 3.13 (probado; se espera que valga $\ge 3.11$) |
| sympy / mpmath | 1.14.0 / 1.3.0 |
| numpy / scipy | 2.5.4 / 1.18.1 |
| matplotlib | 3.11.2 |
| Singular (opcional) | 4.3.2 (solo para `08_contraste_singular.py`; si no está, se omite con código 0) |

## Contenido

### `simbolico/` — cálculo exacto

Cada script imprime resultados legibles, usa comprobaciones explícitas y termina con código 0 solo si todas pasan.
`comun.py` y `explosion.py` son bibliotecas auxiliares (no se ejecutan).

| Script | Qué hace |
|---|---|
| `01_nucleo_cuspide.py` | $\ker(x\mapsto t^2,\,y\mapsto t^3)=(y^2-x^3)$ por eliminación (Gröbner lex) sobre $\mathbb Q$, $\mathbb F_2,\mathbb F_3,\mathbb F_5,\mathbb F_7,\mathbb F_{101}$ |
| `02_potencias_maximal.py` | generadores de $I^n$, $\dim k[x,y]/I^n=n(n+1)/2$, $\dim B/\mathfrak m^n=2n-1$, $n=1..8$ |
| `03_funcion_hilbert.py` | función de Hilbert de $\mathrm{gr}_\mathfrak m B$: cúspide, nodo, $A_n$ ($n=1..5$) y contrastes |
| `04_cartas_blowup.py` | cartas $y=xt$, $x=ys$; transformada total y estricta; $\pi^*C=\tilde C+mE$ (cúspide, nodo, tacnodo) |
| `05_resolucion_cuspide.py` | tres blow-ups sucesivos, tabla de multiplicidades, $\pi^*C=\tilde C+2E_1+3E_2+6E_3$, SNC por Gröbner |
| `06_delta.py` | $\delta$ por la fórmula recursiva y por $\dim_{\mathbb Q}\bar{\mathcal O}/\mathcal O$ directo |
| `07_paraguas_whitney.py` | paraguas $x^2=y^2z$: lugar singular, cartas, lisura, punto pinzado |
| `08_contraste_singular.py` | contraste opcional con Singular (a, b) |

### `figuras/` — figuras (PDF vectorial)

Salida en `figuras/figuras_salida/`: `cuspide_nodo.pdf`, `divisor_proyeccion.pdf`,
`resolucion_cuspide_tres_pasos.pdf` (scripts `fig_*.py`; estilo común en `estilo.py`). Los retratos de fase
(`dinamica_antes.pdf`, `dinamica_despues.pdf`) los produce `dinamica/02_blowup_simulacion.py` en la misma carpeta.
Los PDF no llevan fecha de creación, así que las salidas son reproducibles byte a byte con las mismas versiones.

### `dinamica/` — sistema dinámico del último capítulo de la monografía

| Script | Qué hace | Tipo |
|---|---|---|
| `00_verificacion_simbolica.py` | identidades exactas: tabla de equilibrios del eje, Prop. de equilibrios interiores, cartas $W_1,W_2$ del infinito, los dos blow-ups y la tabla de valores propios, pliegues para $Q=2$ | exacto |
| `01_equilibrios_clasificacion.py` | equilibrios del primer cuadrante (raíces de $h$ + contraste `fsolve`), jacobiano, clasificación, no hiperbólicos marcados; escribe `dinamica/salida/tabla_equilibrios.tex` | numérico |
| `02_blowup_simulacion.py` | blow-up direccional del origen (cartas $x=r,y=rv$ y $y=r,x=ru$), `solve_ivp`, validación cruzada, retratos de fase antes y después | numérico |

`sistema.py` define el campo, el jacobiano, las cartas y **todas las tolerancias y la semilla**:
`solve_ivp` con `DOP853`, `rtol=1e-10`, `atol=1e-12`; semilla `numpy.random.default_rng(20261027)`;
no hiperbólico si $\min|\mathrm{Re}\,\lambda|\le10^{-4}$.

> Todo lo de `dinamica/01` y `02` es **evidencia numérica, no demostración**: el punto flotante no decide igualdades
> (p. ej. $\det J=0$) ni raíces múltiples.
