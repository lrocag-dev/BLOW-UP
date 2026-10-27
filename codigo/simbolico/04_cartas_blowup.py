"""d) Cartas del blow-up del origen de A^2: sustituciones y = x t  y  x = y s.
Transformada total y estricta de la cúspide, el nodo y el tacnodo; se comprueba
pi^* C = C~ + m E con m = multiplicidad, C~ ∩ E = raíces de f_m(1,t) y
(C~·E) = m.
"""
from sympy import symbols, Poly, expand, cancel, groebner, roots, S

from comun import encabezado, comprobar
from explosion import potencia_de, es_unidad

x, y, t, s = symbols("x y t s")

CURVAS = [
    ("cúspide  y^2 = x^3", y**2 - x**3, 2),
    ("nodo  y^2 = x^2 (x+1)", y**2 - x**2 * (x + 1), 2),
    ("tacnodo  y^2 = x^4", y**2 - x**4, 2),
]


def partes(f):
    P = Poly(f, x, y)
    m = min(sum(e) for e in P.monoms())
    fm = sum(c * x**e[0] * y**e[1] for e, c in P.terms() if sum(e) == m)
    return m, fm


def puntos_singulares(g, v1, v2):
    G = groebner([g, g.diff(v1), g.diff(v2)], v1, v2, order="lex")
    return None if es_unidad(G) else G.exprs


if __name__ == "__main__":
    encabezado("d) Cartas del blow-up del origen de A^2")
    for nombre, f, m_esp in CURVAS:
        print(f"\n--- {nombre} ---")
        m, fm = partes(f)
        comprobar(m == m_esp, f"multiplicidad m = {m} (esperada {m_esp}); forma inicial f_m = {fm}")

        # Carta 1: y = x t   (E = {x = 0})
        tot1 = expand(f.subs(y, x * t))
        k1 = potencia_de(tot1, x)
        est1 = expand(cancel(tot1 / x**k1))
        print(f"  carta (x,t), y = x t:  total = {tot1} = x^{k1} * ({est1})")
        comprobar(k1 == m and expand(tot1 - x**m * est1) == 0, "pi^*C = C~ + m E en la carta y = x t")
        comprobar(est1.subs(x, 0) != 0, "x no divide a la transformada estricta (E no es componente de C~)")
        # fórmula f~ = f_m(1,t) + x f_{m+1}(1,t) + ...
        P = Poly(f, x, y)
        formula = sum(c * x ** (sum(e) - m) * t ** e[1] for e, c in P.terms())
        comprobar(expand(est1 - formula) == 0, "f~(x,t) = f_m(1,t) + x f_{m+1}(1,t) + ...")
        # C~ ∩ E : raíces de f_m(1,t)
        corte = est1.subs(x, 0)
        comprobar(expand(corte - fm.subs({x: 1, y: t})) == 0, "C~ ∩ E = {f_m(1,t) = 0}")
        raices = roots(Poly(corte, t))
        grado_inf = m - Poly(corte, t).degree()  # multiplicidad de t = ∞ (dirección x = 0)
        total = sum(raices.values()) + grado_inf
        print(f"  C~ ∩ E en la carta 1: t ∈ {dict(raices)}; multiplicidad en t=∞: {grado_inf}")
        comprobar(total == m, f"(C~·E) = suma de multiplicidades = {total} = m")
        sing1 = puntos_singulares(est1, x, t)
        print(f"  puntos singulares de C~ en la carta 1: {sing1 if sing1 else 'ninguno (C~ lisa)'}")

        # Carta 2: x = y s   (E = {y = 0})
        tot2 = expand(f.subs(x, y * s))
        k2 = potencia_de(tot2, y)
        est2 = expand(cancel(tot2 / y**k2))
        print(f"  carta (y,s), x = y s:  total = {tot2} = y^{k2} * ({est2})")
        comprobar(k2 == m and expand(tot2 - y**m * est2) == 0, "pi^*C = C~ + m E en la carta x = y s")
        print(f"     C~ ∩ E en la carta 2: {est2.subs(y, 0)} = 0")
        comprobar(expand(est2.subs(y, 0) - fm.subs({x: s, y: 1})) == 0,
                  "en la carta x = y s: C~ ∩ E = {f_m(s,1) = 0}")
        if "tacnodo" in nombre:
            comprobar(sorted(map(str, sing1)) == ["t", "x"],
                      "C~ tiene un único punto singular, (x,t) = (0,0): un nodo (t^2 = x^2)")
            comprobar(Poly(est1, x, t).total_degree() == 2 and
                      min(sum(e) for e in Poly(est1, x, t).monoms()) == 2,
                      "mult_(0,0)(C~) = 2: segundo blow-up necesario")
        else:
            comprobar(sing1 is None, "C~ es lisa en la carta 1")
    print("\nOK: pi^*C = C~ + m E con m = 2 en los tres ejemplos.")
