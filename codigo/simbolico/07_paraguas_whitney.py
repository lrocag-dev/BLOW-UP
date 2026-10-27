"""g) Paraguas de Whitney X: x^2 = y^2 z.

  * lugar singular de X (eje z);
  * blow-up de la recta L = V(x,y) (centro liso en A^3, c = 2): dos cartas,
    transformada total y estricta;
  * la transformada estricta es LISA en ambas cartas: no queda singularidad;
  * lo que sí ocurre en z = 0 es el punto pinzado: el morfismo X~ -> X (la
    normalización) deja de ser inmersión exactamente en un punto;
  * el divisor total pi^*X = X~ + 2 E es de cruces normales simples.
"""
from sympy import symbols, groebner, Matrix, expand, cancel, factor, Poly, solve, diff

from comun import encabezado, comprobar
from explosion import potencia_de, es_unidad

x, y, z, s, t = symbols("x y z s t")
f = x**2 - y**2 * z


def reducido_a_ideal(G, base):
    """True si todos los elementos de G pertenecen al ideal (base) (por reducción)."""
    return all(groebner(base, x, y, z, order="grevlex").reduce(g)[1] == 0 for g in G)


if __name__ == "__main__":
    encabezado("g) Paraguas de Whitney  x^2 = y^2 z")

    # --- lugar singular
    J = groebner([f, f.diff(x), f.diff(y), f.diff(z)], x, y, z, order="grevlex")
    print(f"  Ideal jacobiano (f, f_x, f_y, f_z): base de Gröbner = {J.exprs}")
    # J ⊂ (x,y)  y  (x,y)^2 ⊂ J  => rad(J) = (x,y): el eje z
    comprobar(all(groebner([x, y], x, y, z).reduce(g)[1] == 0 for g in J.exprs),
              "J ⊂ (x,y)")
    comprobar(J.reduce(x)[1] == 0 and J.reduce(y**2)[1] == 0 and J.reduce(y * z)[1] == 0,
              "x, y^2, y z ∈ J; luego rad(J) ∋ x, y: Sing X = V(x,y) = eje z")

    # --- blow-up de L = V(x,y)
    # Carta A: x = s y  (variables s, y, z; E = {y = 0})
    totA = expand(f.subs(x, s * y))
    kA = potencia_de(totA, y)
    estA = expand(cancel(totA / y**kA))
    print(f"\n  Carta A (x = s y):  f = {factor(totA)}   ->  X~ : {estA} = 0,   pi^*X = X~ + {kA} E")
    comprobar(kA == 2 and expand(estA - (s**2 - z)) == 0, "carta A: transformada estricta  s^2 = z")
    # Carta B: y = t x  (variables t, x, z; E = {x = 0})
    totB = expand(f.subs(y, t * x))
    kB = potencia_de(totB, x)
    estB = expand(cancel(totB / x**kB))
    print(f"  Carta B (y = t x):  f = {factor(totB)}   ->  X~ : {estB} = 0,   pi^*X = X~ + {kB} E")
    comprobar(kB == 2 and expand(estB - (1 - t**2 * z)) == 0, "carta B: transformada estricta  1 = t^2 z")

    # --- X~ liso en ambas cartas
    GA = groebner([estA, estA.diff(s), estA.diff(y), estA.diff(z)], s, y, z)
    GB = groebner([estB, estB.diff(t), estB.diff(x), estB.diff(z)], t, x, z)
    comprobar(es_unidad(GA), "carta A: X~ es lisa (gradiente = (2s, 0, -1) nunca se anula)")
    comprobar(es_unidad(GB), "carta B: X~ es lisa (el ideal (g, grad g) es (1))")
    print("  => la transformada estricta NO tiene singularidades: no queda singularidad que resolver.")

    # --- divisor total SNC: X~ y E transversales (gradientes independientes)
    JA = Matrix([[estA.diff(s), estA.diff(y), estA.diff(z)], [0, 1, 0]])
    menoresA = [JA[:, [i, j]].det() for i in range(3) for j in range(i + 1, 3)]
    GAE = groebner([estA, y, *menoresA], s, y, z)
    JB = Matrix([[estB.diff(t), estB.diff(x), estB.diff(z)], [0, 1, 0]])
    menoresB = [JB[:, [i, j]].det() for i in range(3) for j in range(i + 1, 3)]
    GBE = groebner([estB, x, *menoresB], t, x, z)
    comprobar(es_unidad(GAE) and es_unidad(GBE), "X~ y E se cortan transversalmente: pi^*X = X~ + 2E es SNC")

    # --- E ∩ X~ y el punto pinzado
    print("\n  E ∩ X~ en la carta A: y = 0, z = s^2  (parábola); pi(s,0,s^2) = (0,0,s^2) en el eje z.")
    print("  El morfismo X~ -> L es 2:1 (s = ±sqrt(z)) y se ramifica en s = 0.")
    comprobar(diff(s**2, s).subs(s, 0) == 0 and Poly(s**2 - symbols("z0"), s).degree() == 2,
              "E∩X~ -> L es de grado 2, ramificado solo en z = 0")
    # La normalización pi: (s,y) |-> (s y, y, s^2) (carta A, superficie X~ parametrizada por (s,y)).
    P = Matrix([s * y, y, s**2])
    Jpi = P.jacobian([s, y])
    menores = [Jpi[[i, j], :].det() for i in range(3) for j in range(i + 1, 3)]
    print(f"  Jacobiana de pi|X~ en la carta A: {Jpi.tolist()}, menores 2x2: {[factor(m) for m in menores]}")
    Gp = groebner(menores, s, y, order="lex")
    comprobar(sorted(map(str, Gp.exprs)) == ["s", "y"],
              "pi|X~ deja de ser inmersión exactamente en (s,y) = (0,0): el punto pinzado, imagen (0,0,0)")
    # En la carta B, pi(t,x) = (x, t x, 1/t^2) es inmersión en todo punto (t != 0)
    Q = Matrix([x, t * x, 1 / t**2])
    JB2 = Q.jacobian([t, x])
    menB = [JB2[[i, j], :].det() for i in range(3) for j in range(i + 1, 3)]
    print(f"  Carta B: menores 2x2 de la jacobiana de pi|X~: {[cancel(m) for m in menB]}")
    comprobar(cancel(menB[1] - 2 / t**3) == 0,
              "carta B: un menor es 2/t^3 != 0 (t != 0 en X~): pi|X~ es inmersión en toda la carta B")
    print("\nOK: Sing(X) = eje z; X~ lisa; punto pinzado = origen; pi^*X = X~ + 2E es SNC.")
