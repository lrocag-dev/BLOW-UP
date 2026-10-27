"""Blow-up de un punto del plano afín, por cartas, con aritmética exacta.

Convenciones (como en la monografía, Cap. 'Curvas planas'):
  carta 1:  b - b0 = (a - a0)*v     ->  a = a0 + u,    b = b0 + u*v
  carta 2:  a - a0 = (b - b0)*v     ->  a = a0 + u*v,  b = b0 + u
En ambas, E = {u = 0}.  Un divisor se representa como lista de componentes
(nombre, polinomio irreducible, coeficiente) sobre las variables de la carta.
"""
from sympy import Poly, symbols, expand, cancel, groebner, Matrix, factor_list, S


def multiplicidad(g, variables, punto):
    """Multiplicidad (orden de anulación) de g en el punto dado."""
    a, b = variables
    u, v = symbols("_u _v")
    h = expand(g.subs({a: punto[0] + u, b: punto[1] + v}, simultaneous=True))
    if h == 0:
        return float("inf")
    return min(sum(e) for e in Poly(h, u, v).monoms())


def potencia_de(g, var):
    """Mayor k con var**k | g (g polinomio no nulo)."""
    P = Poly(expand(g), var)
    return min(e[0] for e in P.monoms())


def explotar(F, componentes, variables, punto, carta, nuevas=("u", "v"), nombre_E="E"):
    """Hace el blow-up de `punto` y devuelve los datos en la carta pedida.

    F: transformada total (polinomio) hasta ahora.
    Devuelve (F_nuevo, componentes_nuevas, variables_nuevas, mult_en_punto, coef_E).
    `mult_en_punto` es un dict nombre -> multiplicidad en el centro.
    """
    a, b = variables
    u, v = symbols(nuevas[0] + " " + nuevas[1])
    if carta == 1:
        sust = {a: punto[0] + u, b: punto[1] + u * v}
    else:
        sust = {a: punto[0] + u * v, b: punto[1] + u}
    mult = {}
    nuevas_comp = []
    coef_E = 0
    for nombre, g, c in componentes:
        m = multiplicidad(g, variables, punto)
        mult[nombre] = m
        coef_E += c * m
        h = expand(g.subs(sust, simultaneous=True))
        k = potencia_de(h, u)
        assert k == m, "la potencia de u extraída debe ser la multiplicidad"
        estricta = cancel(h / u**k)
        if estricta.free_symbols:  # si es constante, la componente no aparece en esta carta
            nuevas_comp.append((nombre, expand(estricta), c))
    nuevas_comp.append((nombre_E, u, coef_E))
    F_nuevo = expand(F.subs(sust, simultaneous=True))
    return F_nuevo, nuevas_comp, (u, v), mult, coef_E


def es_constante_no_nula(e):
    e = cancel(e)
    return (not e.free_symbols) and e != 0


def comprobar_factorizacion(F, componentes):
    """F = c * prod g_i^{e_i} con c constante no nula (la carta es A^2)."""
    prod = S(1)
    for _, g, c in componentes:
        prod *= g**c
    return es_constante_no_nula(cancel(F / prod))


def es_unidad(G):
    return len(G.exprs) == 1 and G.exprs[0] == 1


def es_snc(componentes, variables):
    """Cruces normales simples en la carta (exacto, vía bases de Gröbner):

      1) cada componente es lisa:           (g, g_a, g_b) = (1)
      2) dos componentes se cortan transversalmente:  (g_i, g_j, det Jac) = (1)
      3) no hay puntos triples:             (g_i, g_j, g_k) = (1)
    Devuelve (booleano, lista de mensajes).
    """
    a, b = variables
    polis = [(n, g) for n, g, _ in componentes]
    ok, msgs = True, []
    for n, g in polis:
        G = groebner([g, g.diff(a), g.diff(b)], a, b, order="lex")
        liso = es_unidad(G)
        ok &= liso
        msgs.append(f"{n} lisa: {liso}")
    for i in range(len(polis)):
        for j in range(i + 1, len(polis)):
            (ni, gi), (nj, gj) = polis[i], polis[j]
            det = Matrix([[gi.diff(a), gi.diff(b)], [gj.diff(a), gj.diff(b)]]).det()
            G = groebner([gi, gj, det], a, b, order="lex")
            trans = es_unidad(G)
            ok &= trans
            msgs.append(f"{ni}∩{nj} transversal (o vacía): {trans}")
    for i in range(len(polis)):
        for j in range(i + 1, len(polis)):
            for k in range(j + 1, len(polis)):
                G = groebner([polis[i][1], polis[j][1], polis[k][1]], a, b, order="lex")
                libre = es_unidad(G)
                ok &= libre
                msgs.append(f"{polis[i][0]}∩{polis[j][0]}∩{polis[k][0]} vacía: {libre}")
    return ok, msgs


def numero_interseccion_puntos(gi, gj, variables, extra=()):
    """dim_Q Q[a,b]/(gi,gj,extra...) (número de puntos con multiplicidad)."""
    from comun import dim_cociente
    return dim_cociente([gi, gj, *extra], variables=variables, orden="grevlex")
