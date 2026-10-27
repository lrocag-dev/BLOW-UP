"""c) Función de Hilbert de gr_m(B), B = k[x,y]/(f), para la cúspide, el nodo
y^2 = x^2 (x+1), las singularidades A_n: y^2 = x^(n+1), n = 1..5, y otros
ejemplos de contraste (punto triple ordinario, y^3 = x^4, y^3 = x^5).

H(j) = dim m^j/m^(j+1) = dim B/m^(j+1) - dim B/m^j, con
dim B/m^n = dim k[x,y]/(f, m^n)  (m^n es m-primario, así que el anillo local
y el anillo afín dan la misma dimensión).
Se compara con k[x,y]/(f_m) (la parte inicial) y con la fórmula min(j+1, m).
"""
from sympy import Poly

from comun import (encabezado, comprobar, dim_cociente,
                   generadores_potencia_maximal, x, y)

JMAX = 9  # H(0..JMAX)


def multiplicidad_y_forma_inicial(f):
    P = Poly(f, x, y)
    m = min(sum(e) for e in P.monoms())
    fm = sum(c * x**e[0] * y**e[1] for e, c in P.terms() if sum(e) == m)
    return m, fm


def hilbert_gr(f, jmax=JMAX, modulus=None):
    dims = [0] + [dim_cociente([f] + generadores_potencia_maximal(n), modulus=modulus)
                  for n in range(1, jmax + 2)]
    return [dims[j + 1] - dims[j] for j in range(jmax + 1)]


def hilbert_forma_inicial(fm, jmax=JMAX):
    """Hilbert de k[x,y]/(f_m): dim_j = dim k[x,y]_j - dim (f_m)_j."""
    dims = [0] + [dim_cociente([fm] + generadores_potencia_maximal(n)) for n in range(1, jmax + 2)]
    return [dims[j + 1] - dims[j] for j in range(jmax + 1)]


EJEMPLOS = [
    ("cúspide y^2 = x^3", y**2 - x**3),
    ("nodo y^2 = x^2 (x+1)", y**2 - x**2 * (x + 1)),
    ("A_1: y^2 = x^2", y**2 - x**2),
    ("A_2: y^2 = x^3", y**2 - x**3),
    ("A_3: y^2 = x^4", y**2 - x**4),
    ("A_4: y^2 = x^5", y**2 - x**5),
    ("A_5: y^2 = x^6", y**2 - x**6),
    ("punto triple ordinario xy(x+y)", x * y * (x + y)),
    ("y^3 = x^4", y**3 - x**4),
    ("y^3 = x^5", y**3 - x**5),
]

if __name__ == "__main__":
    encabezado("c) Función de Hilbert de gr_m(B)")
    print(f"{'curva':<34} {'m':>2}  H(0..{JMAX})")
    print("-" * 72)
    for nombre, f in EJEMPLOS:
        m, fm = multiplicidad_y_forma_inicial(f)
        H = hilbert_gr(f)
        esperado = [min(j + 1, m) for j in range(JMAX + 1)]
        print(f"{nombre:<34} {m:>2}  {H}")
        assert H == esperado, f"{nombre}: H = {H} != min(j+1,m) = {esperado}"
        assert H == hilbert_forma_inicial(fm), f"{nombre}: H(gr) != H(k[x,y]/(f_m))"
    # Cúspide: esperado 1,2,2,2,...
    H = hilbert_gr(y**2 - x**3)
    comprobar(H == [1] + [2] * JMAX, "cúspide: H = 1,2,2,2,...")
    H = hilbert_gr(y**2 - x**2 * (x + 1))
    comprobar(H == [1] + [2] * JMAX, "nodo: H = 1,2,2,2,... (m = 2)")
    for n in range(1, 6):
        comprobar(hilbert_gr(y**2 - x ** (n + 1)) == [1] + [2] * JMAX, f"A_{n}: H = 1,2,2,2,...")
    comprobar(True, "H(gr_m B)(j) = min(j+1, m) y coincide con k[x,y]/(f_m) en todos los ejemplos")
