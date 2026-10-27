"""e) Resolución embebida de la cúspide y^2 = x^3 con tres blow-ups sucesivos.

Se calcula, carta por carta y en aritmética exacta:
  * la tabla de multiplicidades (centro, mult de C~, mult del divisor total);
  * pi^*C = C~ + 2 E1 + 3 E2 + 6 E3  (E_i = transformadas estrictas);
  * que el divisor total es de cruces normales simples (SNC) en las cartas
    que recubren la superficie final;
  * autointersecciones E1^2=-3, E2^2=-2, E3^2=-1 y pi^*C . E_i = 0.

Cobertura de la superficie final (4 cartas):
  S1B = carta x = y s del 1.er blow-up (no se vuelve a tocar);
  S2B = carta 'b' del 2.º blow-up; S3A, S3B = las dos cartas del 3.º.
"""
from sympy import symbols, expand

from comun import encabezado, comprobar
from explosion import explotar, comprobar_factorizacion, es_snc, multiplicidad, numero_interseccion_puntos

x, y = symbols("x y")
f = y**2 - x**3


def mostrar(F, comps):
    return " * ".join(f"({g})^{c}" if c != 1 else f"({g})" for _, g, c in comps)


if __name__ == "__main__":
    encabezado("e) Resolución embebida de la cúspide y^2 = x^3 (tres blow-ups)")
    filas = []  # (paso, centro, mult C~, mult divisor total ponderado, coef nuevo E)

    # ---------------- Paso 1: blow-up de (0,0) en A^2 ----------------
    F0, comps0 = f, [("C~", f, 1)]
    F1a, c1a, v1a, mult1, coefE1 = explotar(F0, comps0, (x, y), (0, 0), 1, ("x", "t"), "E1")
    F1b, c1b, v1b, _, _ = explotar(F0, comps0, (x, y), (0, 0), 2, ("y", "s"), "E1")
    filas.append(("1", "(0,0) ∈ A^2", mult1["C~"], "-", "-", coefE1))
    print(f"\nPaso 1 | carta 1 (y = x t): F = {mostrar(F1a, c1a)}")
    print(f"       | carta 2 (x = y s): F = {mostrar(F1b, c1b)}")
    comprobar(coefE1 == 2, "coeficiente de E1 = mult_0(C) = 2")
    comprobar(comprobar_factorizacion(F1a, c1a) and comprobar_factorizacion(F1b, c1b),
              "pi_1^*C = C~ + 2 E1 en ambas cartas")
    t_, x_ = v1a[1], v1a[0]  # (x, t)
    C1, Ea1 = c1a[0][1], c1a[1][1]
    ok1a, _ = es_snc(c1a, v1a)
    comprobar(not ok1a, "carta (x,t): C~ (t^2 = x) y E1 (x = 0) son tangentes: el divisor NO es SNC")
    ok, msg = es_snc(c1b, v1b)
    comprobar(ok, "S1B (carta x = y s): divisor SNC (solo E1; C~ no corta a E1 aquí)")

    # ---------------- Paso 2: blow-up de q1 = C~ ∩ E1 ----------------
    q1 = (0, 0)  # en la carta (x,t)
    comprobar(C1.subs({v1a[0]: 0, v1a[1]: 0}) == 0 and Ea1.subs({v1a[0]: 0, v1a[1]: 0}) == 0,
              "q1 = (x,t) = (0,0) está en C~ y en E1 (C~ es lisa y tangente a E1: t^2 = x)")
    # La dirección tangente común es x = 0 (E1), que está en la carta 2: x = t w.
    F2a, c2a, v2a, mult2, coefE2 = explotar(F1a, c1a, v1a, q1, 2, ("t", "w"), "E2")
    F2b, c2b, v2b, _, _ = explotar(F1a, c1a, v1a, q1, 1, ("x", "tau"), "E2")
    filas.append(("2", "q1 = C~∩E1", mult2["C~"], mult2["E1"], "-", coefE2))
    print(f"\nPaso 2 | carta (t,w), x = t w:   F = {mostrar(F2a, c2a)}")
    print(f"       | carta (x,tau), t = x tau: F = {mostrar(F2b, c2b)}")
    comprobar(mult2 == {"C~": 1, "E1": 1}, f"multiplicidades en q1: {mult2}")
    comprobar(coefE2 == 3, "coeficiente de E2 = 1*1 + 2*1 = 3")
    comprobar(comprobar_factorizacion(F2a, c2a) and comprobar_factorizacion(F2b, c2b),
              "factorización del transformado total en ambas cartas")
    ok, _ = es_snc(c2b, v2b)
    comprobar(ok, "S2B (carta t = x tau): divisor SNC")
    comprobar(len(c2a) == 3, "carta (t,w): pasan C~, E1 y E2 por el origen (tres ramas): NO es SNC todavía")
    ok2a, _ = es_snc(c2a, v2a)
    comprobar(not ok2a, "carta (t,w): el divisor NO es SNC (tres curvas lisas concurrentes)")

    # ---------------- Paso 3: blow-up de q2 (punto triple) ----------------
    q2 = (0, 0)  # en la carta (t,w)
    F3a, c3a, v3a, mult3, coefE3 = explotar(F2a, c2a, v2a, q2, 1, ("t", "s"), "E3")
    F3b, c3b, v3b, _, _ = explotar(F2a, c2a, v2a, q2, 2, ("w", "tau"), "E3")
    filas.append(("3", "q2 = C~∩E1∩E2", mult3["C~"], mult3["E1"], mult3["E2"], coefE3))
    print(f"\nPaso 3 | carta (t,s), w = t s:   F = {mostrar(F3a, c3a)}")
    print(f"       | carta (w,tau), t = w tau: F = {mostrar(F3b, c3b)}")
    comprobar(mult3 == {"C~": 1, "E1": 1, "E2": 1}, f"multiplicidades en q2: {mult3}")
    comprobar(coefE3 == 6, "coeficiente de E3 = 1*1 + 2*1 + 3*1 = 6")
    comprobar(comprobar_factorizacion(F3a, c3a) and comprobar_factorizacion(F3b, c3b),
              "factorización del transformado total en ambas cartas")
    coef = {n: c for n, _, c in c3a}
    coef.update({n: c for n, _, c in c3b})
    comprobar(coef == {"C~": 1, "E1": 2, "E2": 3, "E3": 6},
              f"pi^*C = C~ + 2 E1 + 3 E2 + 6 E3   (coeficientes calculados: {coef})")

    # ---------------- SNC final ----------------
    for nombre, comps, vs in [("S1B", c1b, v1b), ("S2B", c2b, v2b), ("S3A", c3a, v3a), ("S3B", c3b, v3b)]:
        ok, msgs = es_snc(comps, vs)
        comprobar(ok, f"{nombre}: divisor SNC  [componentes: {[n for n,_,_ in comps]}]")

    # ---------------- Intersecciones ----------------
    # Autointersecciones: E_i^2 = -1 - #(centros posteriores que pertenecen a E_i~)
    centros = [("q1", {"C~": mult2["C~"], "E1": mult2["E1"]}), ("q2", mult3)]
    autoint = {"E1": -1, "E2": -1, "E3": -1}
    for _, m in centros:
        for n in ("E1", "E2"):
            if m.get(n, 0) >= 1:
                autoint[n] -= 1
    comprobar(autoint == {"E1": -3, "E2": -2, "E3": -1}, f"autointersecciones {autoint}")

    # Intersecciones entre componentes distintas (puntos con multiplicidad):
    # en S3A todos los puntos; en S3B solo los de tau = 0 (el resto ya está en S3A).
    def inter(comps, vs, i, j, extra=()):
        d = {n: g for n, g, _ in comps}
        return numero_interseccion_puntos(d[i], d[j], vs, extra)

    M = {}
    for (i, j) in [("E3", "E1"), ("E3", "C~"), ("E1", "C~")]:
        M[(i, j)] = inter(c3a, v3a, i, j)
    tau = v3b[1]
    M[("E3", "E2")] = inter(c3b, v3b, "E3", "E2", (tau,))
    M[("E2", "C~")] = inter(c3b, v3b, "E2", "C~", (tau,))
    M[("E1", "E2")] = 0  # E1 y E2 no aparecen juntas en S3A ni en S3B; S1B y S2B tienen una sola E
    print("\n  Intersecciones (puntos, transversales):", M)
    comprobar(M[("E3", "E1")] == M[("E3", "E2")] == M[("E3", "C~")] == 1
              and M[("E1", "C~")] == M[("E2", "C~")] == M[("E1", "E2")] == 0,
              "grafo dual: E3 corta a C~, E1 y E2 (1 punto cada una); no hay otras intersecciones")
    coefs = {"C~": 1, "E1": 2, "E2": 3, "E3": 6}
    def par(i, j):
        if i == j:
            return autoint[i]
        return M.get((i, j), M.get((j, i), 0))
    for Ei in ("E1", "E2", "E3"):
        valor = sum(c * par(n, Ei) if n != "C~" else c * M.get((Ei, "C~"), M.get(("C~", Ei), 0))
                    for n, c in coefs.items())
        comprobar(valor == 0, f"pi^*C . {Ei} = 0  (sumando coef * (componente . {Ei}))")

    # ---------------- Tabla de multiplicidades ----------------
    print("\nTabla de multiplicidades de la sucesión de blow-ups:")
    print(f"{'paso':>4} | {'centro':<15} | {'m(C~)':>5} | {'m(E1~)':>6} | {'m(E2~)':>6} | coef. nuevo E = sum coef*mult")
    print("-" * 78)
    for paso, centro, mC, m1, m2, ce in filas:
        print(f"{paso:>4} | {centro:<15} | {mC:>5} | {m1:>6} | {m2:>6} | {ce}")
    comprobar([r[2] for r in filas] == [2, 1, 1] and [r[5] for r in filas] == [2, 3, 6],
              "multiplicidades de C~: (2,1,1);  coeficientes de E1,E2,E3: (2,3,6)")
    delta = sum(r[2] * (r[2] - 1) // 2 for r in filas)
    comprobar(delta == 1, f"delta = sum m_i(m_i-1)/2 = {delta}")
    print("\nOK: resolución embebida de la cúspide verificada.")
