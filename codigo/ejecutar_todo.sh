#!/usr/bin/env bash
# Reproduce TODOS los resultados del apéndice computacional en un entorno limpio.
# Uso (desde cualquier directorio):   bash codigo/ejecutar_todo.sh
# Crea codigo/.venv, instala las versiones fijadas de requirements.txt, ejecuta todos
# los scripts (cada uno termina con código 0 solo si pasan todos sus assert) y deja el
# registro completo en codigo/registro_ejecucion.txt.
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -x .venv/bin/python ]; then
  python3 -m venv .venv
fi
.venv/bin/python -m pip install --quiet --upgrade pip
.venv/bin/python -m pip install --quiet -r requirements.txt

SCRIPTS=(
  simbolico/01_nucleo_cuspide.py
  simbolico/02_potencias_maximal.py
  simbolico/03_funcion_hilbert.py
  simbolico/04_cartas_blowup.py
  simbolico/05_resolucion_cuspide.py
  simbolico/06_delta.py
  simbolico/07_paraguas_whitney.py
  simbolico/08_contraste_singular.py   # opcional: se omite si no hay Singular
  figuras/fig_cuspide_nodo.py
  figuras/fig_divisor_proyeccion.py
  figuras/fig_resolucion_tres_pasos.py
  dinamica/00_verificacion_simbolica.py
  dinamica/01_equilibrios_clasificacion.py
  dinamica/02_blowup_simulacion.py
)

REGISTRO=registro_ejecucion.txt
{
  echo "# Registro de ejecución (generado por codigo/ejecutar_todo.sh)"
  echo "# $(.venv/bin/python --version); $(.venv/bin/python -m pip freeze | grep -iE '^(sympy|mpmath|numpy|scipy|matplotlib)==' | tr '\n' ' ')"
  for s in "${SCRIPTS[@]}"; do
    echo
    echo "################ $s ################"
    .venv/bin/python "$s"
    echo "[código de salida 0] $s"
  done
  echo
  echo "TODOS LOS SCRIPTS TERMINARON CON CÓDIGO 0"
} 2>&1 | tee "$REGISTRO"
