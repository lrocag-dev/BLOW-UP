"""a) ker(phi: k[x,y] -> k[t], x -> t^2, y -> t^3) = (y^2 - x^3).

Teorema de eliminación: con el orden lexicográfico t > x > y, la intersección
de la base de Gröbner de (x - t^2, y - t^3) con k[x,y] genera el núcleo.
Se hace sobre Q y sobre varios cuerpos finitos F_p.
"""
from sympy import symbols, groebner, expand

from comun import encabezado, comprobar

t, x, y = symbols("t x y")


def nucleo(modulus=None):
    G = groebner([x - t**2, y - t**3], t, x, y, order="lex", modulus=modulus)
    return [g for g in G.exprs if not g.has(t)], G


if __name__ == "__main__":
    encabezado("a) Núcleo de phi: x -> t^2, y -> t^3")
    comprobar(expand((t**3) ** 2 - (t**2) ** 3) == 0, "phi(y^2 - x^3) = t^6 - t^6 = 0")

    esperado = groebner([y**2 - x**3], x, y, order="lex").exprs  # [x^3 - y^2]
    for cuerpo, p in [("Q", None), ("F_2", 2), ("F_3", 3), ("F_5", 5), ("F_7", 7), ("F_101", 101)]:
        elim, G = nucleo(p)
        esp = groebner([y**2 - x**3], x, y, order="lex", modulus=p).exprs
        print(f"  Base de Gröbner lex (t>x>y) sobre {cuerpo}: {G.exprs}")
        print(f"  Parte sin t (ideal de eliminación):      {elim}")
        comprobar(len(elim) == 1, f"el ideal de eliminación sobre {cuerpo} tiene un único generador")
        comprobar(elim == esp, f"ker(phi) = (y^2 - x^3) sobre {cuerpo}")
    print("\nOK: ker(phi) = (y^2 - x^3) en todos los cuerpos probados.")
