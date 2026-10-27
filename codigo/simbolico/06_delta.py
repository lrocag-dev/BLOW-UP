"""f) Invariante delta de curvas planas.

  (i)  Por la fórmula de la monografía  delta = m(m-1)/2 + sum_q delta_q,
       aplicada recursivamente sobre los puntos infinitamente próximos
       (se exige que las direcciones tangentes sean racionales).
  (ii) De forma directa: delta = dim_Q (Obar/O), donde Obar = prod Q[[t_i]]
       (una rama por factor), usando las parametrizaciones exactas de las ramas.
       Se trabaja módulo t^N con N tal que t^N Obar ⊂ O (conductor), de modo que
       delta = sum N_i - dim(imagen de O en Obar/t^N).  Se verifica la
       estabilidad para dos valores de N.
Todo con racionales exactos (sympy / DomainMatrix); nunca flotantes.
"""
from fractions import Fraction

from sympy import symbols, Poly, expand, cancel, QQ, Rational, binomial
from sympy.polys.matrices import DomainMatrix

from comun import encabezado, comprobar

x, y, t = symbols("x y t")


# ------------------------------------------------------------------ (i) fórmula
def forma_inicial(f, v1, v2):
    P = Poly(f, v1, v2)
    m = min(sum(e) for e in P.monoms())
    return m, sum(c * v1**e[0] * v2**e[1] for e, c in P.terms() if sum(e) == m)


def delta_recursivo(f, profundidad=0):
    """delta del germen de f en el origen de (x,y), por la fórmula de la monografía.

    Devuelve (delta, arbol) con arbol = [m, subarbol_1, ...] (multiplicidad en
    el punto y árboles de los puntos de la transformada estricta sobre E).
    """
    assert profundidad < 40, "demasiadas explosiones"
    f = expand(f)
    m, _ = forma_inicial(f, x, y)
    if m == 1:
        return 0, [1]
    b = symbols("b")
    # carta 1: y = x b ; los puntos de E son las raíces de f_m(1, b)
    est1 = expand(cancel(expand(f.subs(y, x * b)) / x**m))
    g1 = Poly(est1.subs(x, 0), b)
    raices = g1.ground_roots()  # raíces racionales con multiplicidad
    n_inf = m - g1.degree()  # multiplicidad de la dirección x = 0 (b = infinito)
    if sum(raices.values()) + n_inf != m:
        raise NotImplementedError("direcciones tangentes no racionales")
    d = m * (m - 1) // 2
    arbol = [m]
    for r in raices:
        germen = est1.subs(b, b + r).subs({x: x, b: y}, simultaneous=True)
        dq, ar = delta_recursivo(germen, profundidad + 1)
        d += dq
        arbol.append(ar)
    if n_inf > 0:
        # carta 2: x = y b ; el punto de E es b = 0 (coordenadas (y, b) -> (x, y))
        est2 = expand(cancel(expand(f.subs(x, y * b)) / y**m))
        germen = est2.subs({y: x, b: y}, simultaneous=True)
        dq, ar = delta_recursivo(germen, profundidad + 1)
        d += dq
        arbol.append(ar)
    return d, arbol


def aplanar(arbol):
    if isinstance(arbol, int):
        return [arbol]
    out = [arbol[0]]
    for h in arbol[1:]:
        out += aplanar(h)
    return out


# ------------------------------------------------------------------ (ii) directo
def a_lista(P, N):
    """Polinomio en t -> lista de N coeficientes racionales (Fraction)."""
    P = Poly(P, t)
    out = []
    for k in range(N):
        c = P.coeff_monomial(t**k)
        out.append(Fraction(int(c.p), int(c.q)))
    return out


def mult_serie(A, B, N):
    C = [Fraction(0)] * N
    for i, a in enumerate(A):
        if a:
            for j in range(N - i):
                C[i + j] += a * B[j]
    return C


def serie_raiz_1_mas_t(N):
    """sqrt(1+t) truncado en t^N, con coeficientes binomial(1/2, k) exactos."""
    return sum(binomial(Rational(1, 2), k) * t**k for k in range(N))


def delta_directo(ramas, N):
    """ramas: lista de pares (X(t), Y(t)) con coeficientes racionales.

    delta = sum N_i - rango de {x^i y^j} en prod Q[t]/t^N.  Como ord_t(x),
    ord_t(y) >= 1 en cada rama, los monomios con i + j >= N son nulos mod t^N,
    así que basta i + j <= N - 1.
    """
    n = len(ramas)
    datos = []
    for X, Y in ramas:
        A, B = a_lista(X, N), a_lista(Y, N)
        unidad = [Fraction(1)] + [Fraction(0)] * (N - 1)
        px, py = [unidad], [unidad]
        for _ in range(N - 1):
            px.append(mult_serie(px[-1], A, N))
            py.append(mult_serie(py[-1], B, N))
        datos.append((px, py))
    filas = []
    for i in range(N):
        for j in range(N - i):
            fila = []
            for px, py in datos:
                fila += mult_serie(px[i], py[j], N)
            filas.append([QQ(c.numerator, c.denominator) for c in fila])
    rango = DomainMatrix(filas, (len(filas), n * N), QQ).rank()
    return n * N - rango


# ------------------------------------------------------------------ ejemplos
def A_n(n):
    f = y**2 - x ** (n + 1)
    if (n + 1) % 2 == 1:  # unirrama
        ramas = [(t**2, t ** (n + 1))]  # x = t^2, y = t^(n+1)
    else:
        k = (n + 1) // 2
        ramas = [(t, t**k), (t, -(t**k))]
    return f, ramas


def ejemplos(N):
    ej = []
    for n in range(1, 6):
        f, ramas = A_n(n)
        ej.append((f"A_{n}: y^2 = x^{n+1}", f, ramas, (n + 1) // 2))
    ej.append(("cúspide y^2 = x^3", y**2 - x**3, [(t**2, t**3)], 1))
    ej.append(("nodo y^2 = x^2 (x+1)", y**2 - x**2 * (x + 1),
               [(t, t * serie_raiz_1_mas_t(N)), (t, -t * serie_raiz_1_mas_t(N))], 1))
    ej.append(("tacnodo y^2 = x^4", y**2 - x**4, [(t, t**2), (t, -(t**2))], 2))
    ej.append(("punto triple ordinario xy(x+y)", x * y * (x + y), [(t, 0 * t), (0 * t, t), (t, -t)], 3))
    ej.append(("y^3 = x^4  (tipo E_6)", y**3 - x**4, [(t**3, t**4)], 3))
    ej.append(("y^3 = x^5  (tipo E_8)", y**3 - x**5, [(t**3, t**5)], 4))
    return ej


if __name__ == "__main__":
    encabezado("f) delta por la fórmula de la monografía y por cálculo directo de dim(Obar/O)")
    print(f"{'curva':<32} | {'mult. (sucesión)':<22} | {'delta fórmula':>13} | {'delta directo':>13} | esperado")
    print("-" * 96)
    N1, N2 = 12, 16
    for nombre, f, ramas, esperado in ejemplos(N2):
        for X, Y in ramas:  # las parametrizaciones parametrizan realmente a la curva (módulo t^N)
            val = Poly(expand(f.subs({x: X, y: Y})), t)
            assert all(val.coeff_monomial(t**k) == 0 for k in range(N2 - 3)), f"{nombre}: parametrización incorrecta"
        d_formula, arbol = delta_recursivo(f)
        d1 = delta_directo(ramas, N1)
        d2 = delta_directo(ramas, N2)
        seq = [m for m in aplanar(arbol) if m >= 2]
        print(f"{nombre:<32} | {str(seq):<22} | {d_formula:>13} | {d1:>13} | {esperado}")
        assert d_formula == esperado, f"{nombre}: delta (fórmula) = {d_formula} != {esperado}"
        assert d1 == d2 == esperado, f"{nombre}: delta directo = {d1}, {d2} != {esperado}"
        assert sum(m * (m - 1) // 2 for m in seq) == d_formula, 'delta = suma de m(m-1)/2'
    comprobar(True, "delta (fórmula recursiva) = delta (dim Obar/O) en todos los ejemplos, estable en N = 12, 16")
    comprobar(True, "A_n: delta = ceil(n/2):  1,1,2,2,3 para n=1..5")
