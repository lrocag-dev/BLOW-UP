"""Contraste opcional con Singular (si está instalado; si no, se omite con código 0).

Repite con Singular los cálculos (a) y (b): el núcleo de phi por eliminación
y las dimensiones dim k[x,y]/m^n = n(n+1)/2 y dim B/m^n = 2n-1 (n = 1..8),
y compara con lo que dan las bases de Gröbner de sympy.
"""
import shutil
import subprocess
import sys

from sympy import symbols, groebner, sympify, expand

from comun import encabezado, comprobar

NMAX = 8
PROGRAMA = f"""
LIB "elim.lib";
ring r=0,(t,x,y),dp;
short=0;
ideal J=eliminate(ideal(x-t^2,y-t^3),t);
"ELIM";
J[1];
size(J);
"FIN";
ring s=0,(x,y),dp;
ideal f=y^2-x^3;
int n;
for(n=1;n<={NMAX};n++){{ "A", vdim(std(maxideal(n))); "B", vdim(std(f+maxideal(n))); }}
quit;
"""

if __name__ == "__main__":
    encabezado("Contraste con Singular (opcional)")
    exe = shutil.which("Singular")
    if exe is None:
        print("  Singular no está instalado: contraste OMITIDO (no es obligatorio).")
        sys.exit(0)
    salida = subprocess.run([exe, "-q"], input=PROGRAMA, capture_output=True, text=True, timeout=120).stdout
    lineas = [l.strip() for l in salida.splitlines() if l.strip() and "redefining" not in l]
    print("  Salida de Singular:", lineas[:7], "...")
    x, y = symbols("x y")
    i = lineas.index("ELIM")
    gen = sympify(lineas[i + 1].replace("^", "**"))
    esperado = groebner([y**2 - x**3], x, y, order="lex").exprs[0]
    comprobar(lineas[i + 2] == "1", "Singular: el ideal de eliminación tiene un único generador")
    comprobar(expand(gen - esperado) == 0 or expand(gen + esperado) == 0,
              f"Singular: eliminación da {gen}  (sympy: {esperado})")
    A = [int(l.split()[1]) for l in lineas if l.startswith("A ")]
    B = [int(l.split()[1]) for l in lineas if l.startswith("B ")]
    comprobar(A == [n * (n + 1) // 2 for n in range(1, NMAX + 1)], f"Singular: dim k[x,y]/m^n = {A}")
    comprobar(B == [2 * n - 1 for n in range(1, NMAX + 1)], f"Singular: dim B/m^n = {B}")
