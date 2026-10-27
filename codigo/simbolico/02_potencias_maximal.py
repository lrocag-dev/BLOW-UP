"""b) Potencias de I = (x,y): generadores, dim k[x,y]/I^n = n(n+1)/2 y
dim B/m^n = 2n-1 para B = k[x,y]/(y^2 - x^3), n = 1..8.
"""
from sympy import groebner

from comun import (encabezado, comprobar, dim_cociente,
                   generadores_potencia_maximal, x, y)

NMAX = 8
CUERPOS = [("Q", None), ("F_2", 2), ("F_3", 3), ("F_101", 101)]

if __name__ == "__main__":
    encabezado("b) Potencias de I = (x,y) y dimensiones de cocientes")
    print(f"{'n':>2} | {'#gen I^n':>8} | {'dim k[x,y]/I^n':>14} | {'n(n+1)/2':>8} | {'dim B/m^n':>9} | {'2n-1':>5}")
    print("-" * 62)
    for n in range(1, NMAX + 1):
        gens = generadores_potencia_maximal(n)
        # I^n está generado por los n+1 monomios x^i y^(n-i)
        esperados = {x**i * y ** (n - i) for i in range(n + 1)}
        assert set(gens) == esperados, "generadores de I^n"
        G = groebner(gens, x, y, order="grevlex")
        assert set(G.exprs) == esperados, "la base de Gröbner de I^n son los mismos monomios"
        for nombre, p in CUERPOS:
            dA = dim_cociente(gens, modulus=p)
            dB = dim_cociente([y**2 - x**3] + gens, modulus=p)
            assert dA == n * (n + 1) // 2, f"dim k[x,y]/I^{n} sobre {nombre}"
            assert dB == 2 * n - 1, f"dim B/m^{n} sobre {nombre}"
            if p is None:
                print(f"{n:>2} | {len(gens):>8} | {dA:>14} | {n*(n+1)//2:>8} | {dB:>9} | {2*n-1:>5}")
    comprobar(True, f"n=1..{NMAX}: #gen(I^n)=n+1, dim k[x,y]/I^n=n(n+1)/2, dim B/m^n=2n-1 "
                    f"(sobre {', '.join(c for c, _ in CUERPOS)})")
