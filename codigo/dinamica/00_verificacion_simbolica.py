"""Verificación simbólica EXACTA del sistema de la monografía (Cap. de dinámica):

    x' = (f(x) - y)/2,   y' = P y g(x,y),
    f(x) = (1 - x^2) x (x + Q),   g = x - (S + R y)(x + Q).

Comprueba con sympy (sin flotantes): la tabla de equilibrios del eje (a), la
Proposición sobre equilibrios interiores (b), las ecuaciones de la carta W_2 del
infinito, los dos blow-ups direccionales de (c) y la tabla de valores propios.
Los parámetros P, Q, R, S son símbolos positivos.
"""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simbolico"))
from comun import encabezado, comprobar  # noqa: E402

x, y, v, z, zeta, u, w, xi, eta = sp.symbols("x y v z zeta u w xi eta")
P, Q, R, S = sp.symbols("P Q R S", positive=True)

f = (1 - x**2) * x * (x + Q)
g = x - (S + R * y) * (x + Q)
X = (f - y) / 2
Y = P * y * g


def jac(Xc, Yc, a, b):
    return sp.Matrix([[sp.diff(Xc, a), sp.diff(Xc, b)], [sp.diff(Yc, a), sp.diff(Yc, b)]])


def explotar_campo(Xa, Xb, a, b, carta, u_, s_):
    """Blow-up direccional. carta 1: b = a s (E={a=0}); carta 2: a = b s (E={b=0}).
    Devuelve el campo (en u, s) dividido por la potencia máxima de u y esa potencia."""
    if carta == 1:
        sust = {b: a * s_}
        Xu = Xa
        Xs = (Xb - s_ * Xa) / a
        base = a
    else:
        sust = {a: b * s_}
        Xu = Xb
        Xs = (Xa - s_ * Xb) / b
        base = b
    Xu = sp.expand(sp.cancel(Xu.subs(sust, simultaneous=True)))
    Xs = sp.expand(sp.cancel(Xs.subs(sust, simultaneous=True)))
    # renombrar la variable radial a u_
    Xu, Xs = Xu.subs(base, u_), Xs.subs(base, u_)
    k = min(min(m[0] for m in sp.Poly(c, u_).monoms()) for c in (Xu, Xs) if c != 0)
    return sp.expand(Xu / u_**k), sp.expand(Xs / u_**k), k


if __name__ == "__main__":
    encabezado("Verificación simbólica del sistema (P,Q,R,S > 0)")

    # ---------- (a) equilibrios en el eje y = 0 ----------
    print("\n(a) Equilibrios del eje y = 0 (jacobiana triangular):")
    J = jac(X, Y, x, y)
    tabla = {
        0: (Q / 2, -P * Q * S),
        1: (-(1 + Q), P * (1 - S * (1 + Q))),
        -1: (1 - Q, P * (S * (1 - Q) - 1)),
        -Q: (Q * (Q**2 - 1) / 2, -P * Q),
    }
    for x0, (l1, l2) in tabla.items():
        comprobar(sp.simplify((X.subs({x: x0, y: 0}))) == 0 and sp.simplify(Y.subs({x: x0, y: 0})) == 0,
                  f"(x0,0) con x0 = {x0} es equilibrio")
        J0 = J.subs({x: x0, y: 0})
        comprobar(J0[1, 0] == 0, "jacobiana triangular superior en y = 0")
        comprobar(sp.expand(J0[0, 0] - l1) == 0 and sp.expand(J0[1, 1] - l2) == 0,
                  f"valores propios en x0 = {x0}: {sp.factor(l1)}, {sp.factor(l2)}")
    # coincidencia de (-1,0) y (-Q,0) en Q = 1, y valores propios 0, -P
    J_sn = J.subs({x: -1, y: 0, Q: 1})
    comprobar(sp.simplify(J_sn[0, 0]) == 0 and sp.simplify(J_sn[1, 1]) == -P,
              "Q = 1: valores propios 0 y -P en (-1,0)")
    # coeficiente cuadrático sobre el eje en Q = 1: x' = f(x)/2, con x = -1 + e
    e = sp.symbols("e")
    sobre_eje = sp.series((f.subs({Q: 1, x: -1 + e}) / 2), e, 0, 3).removeO()
    comprobar(sp.expand(sobre_eje + e**2) == 0 or sp.expand(sobre_eje - (-e**2)) == 0,
              f"Q = 1, sobre el eje: x' = {sobre_eje} + ... (silla-nodo, coef. -1)")

    # ---------- (b) equilibrios interiores ----------
    print("\n(b) Equilibrios fuera del eje:")
    Phi = x / (x + Q) - R * f
    h = sp.expand(R * (x + Q) * f - (1 - S) * x + S * Q)
    comprobar(sp.simplify((x + Q) * (S - Phi) - h) == 0, "h(x) = (x+Q)(S - Phi(x))")
    comprobar(sp.expand(g.subs(y, f) + h) == 0, "g(x, f(x)) = -h(x)")
    comprobar(sp.Poly(h, x).degree() == 5 and sp.Poly(h, x).LC() == -R, "h tiene grado 5 y coeficiente principal -R")
    for punto, valor in [(0, S * Q), (-Q, Q), (1, S * (1 + Q) - 1), (-1, 1 - S * (1 - Q))]:
        comprobar(sp.expand(h.subs(x, punto) - valor) == 0, f"h({punto}) = {valor}")
    # jacobiana sobre la curva y = f(x):  det J = -(P/2) f h'  y  tr J en h = 0
    Jc = J.subs(y, f)
    Jc = Jc.applyfunc(sp.expand)
    det_c, tr_c = sp.expand(Jc.det()), sp.expand(Jc.trace())
    hp = sp.diff(h, x)
    # Identidad módulo h: (det + (P/2) f h') = 0 mod h  y  det - (P/2) f (x+Q) Phi' = 0 mod h
    num1 = sp.expand(det_c + P / 2 * f * hp)
    r1 = sp.rem(sp.Poly(num1, x), sp.Poly(h, x))
    comprobar(r1.is_zero, "det J = -(P/2) f h'  módulo h(x)   (exacto, resto de la división = 0)")
    # Phi' tiene denominador (x+Q)^2: multiplicar
    Phip = sp.together(sp.diff(Phi, x))
    num2 = sp.expand(sp.cancel((det_c - P / 2 * f * (x + Q) * Phip) * (x + Q) ** 2))
    r2 = sp.rem(sp.Poly(num2, x), sp.Poly(h, x))
    comprobar(r2.is_zero, "det J = (P/2) f (x+Q) Phi'  módulo h(x)")
    num3 = sp.expand(tr_c - (sp.diff(f, x) / 2 - P * R * (x + Q) * f))
    r3 = sp.rem(sp.Poly(num3, x), sp.Poly(h, x))
    comprobar(r3.is_zero, "tr J = f'/2 - P R (x+Q) f  módulo h(x)")

    # ---------- (c) infinito: carta W_1 (direcciones (±1,0)) ----------
    print("\n(c) El infinito, carta W_1: (x,y) = (1/z, u/z), campo multiplicado por z^3:")
    sub1 = {x: 1 / z, y: u / z}
    zdot1 = sp.expand(sp.cancel(-(z ** 5) * X.subs(sub1, simultaneous=True)))
    udot1 = sp.expand(sp.cancel(z**4 * (Y - u * X).subs(sub1, simultaneous=True)))
    print("    z' =", zdot1)
    print("    u' =", udot1)
    Pz = sp.Poly(zdot1, z, u)
    Pu = sp.Poly(udot1, z, u)
    comprobar(sp.expand(Pz.coeff_monomial(z) - sp.Rational(1, 2)) == 0
              and sp.expand(Pz.coeff_monomial(z**2) - Q / 2) == 0
              and sp.expand(Pz.coeff_monomial(z**3) + sp.Rational(1, 2)) == 0,
              "z' = z/2 + (Q/2) z^2 - (1/2) z^3 + ...  (coincide con la monografía)")
    comprobar(sp.expand(Pu.coeff_monomial(u) - sp.Rational(1, 2)) == 0
              and sp.expand(Pu.coeff_monomial(u * z) - Q / 2) == 0
              and sp.expand(Pu.coeff_monomial(u**2 * z) + P * R) == 0,
              "u' = u/2 + (Q/2) u z - PR u^2 z + ...  (coincide con la monografía)")
    J1inf = jac(zdot1, udot1, z, u).subs({z: 0, u: 0})
    comprobar(J1inf == sp.eye(2) / 2, "parte lineal en (z,u) = (0,0): Id/2 (nodo hiperbólico, valor propio doble 1/2)")
    comprobar(sp.expand(udot1.subs(z, 0) - u / 2) == 0 and zdot1.subs(z, 0) == 0,
              "sobre z = 0: u' = u/2, luego (0,0) es el único equilibrio del infinito en W_1")

    # ---------- (c) infinito: carta W_2 ----------
    print("\n(c) El infinito, carta W_2: (x,y) = (v/z, 1/z), campo multiplicado por z^3:")
    d = 4
    sub = {x: v / z, y: 1 / z}
    vdot = sp.expand(sp.cancel(z**d * (X - v * Y).subs(sub, simultaneous=True)))
    zdot = sp.expand(sp.cancel(-z ** (d + 1) * Y.subs(sub, simultaneous=True)))
    Rr, Pp, Qq, Ss = R, P, Q, S
    vdot_mono = (-sp.Rational(1, 2) * z**3 + Pp * Rr * v**2 * z + Pp * Qq * Rr * v * z**2
                 - sp.Rational(1, 2) * v**4 - Qq / 2 * v**3 * z + (Pp * (Ss - 1) + sp.Rational(1, 2)) * v**2 * z**2
                 + (Pp * Qq * Ss + Qq / 2) * v * z**3)
    zdot_mono = Pp * Rr * v * z**2 + Pp * Qq * Rr * z**3 + Pp * (Ss - 1) * v * z**3 + Pp * Qq * Ss * z**4
    print("    v' =", vdot)
    print("    z' =", zdot)
    comprobar(sp.expand(vdot - vdot_mono) == 0, "v' coincide con la ecuación de la monografía")
    comprobar(sp.expand(zdot - zdot_mono) == 0, "z' coincide con la ecuación de la monografía")
    J_inf = jac(vdot, zdot, v, z).subs({v: 0, z: 0})
    comprobar(J_inf == sp.zeros(2, 2), "el origen de W_2 es totalmente degenerado (J = 0)")
    # parte homogénea de grado 3 y polinomio característico
    def parte_grado(expr, k):
        Pl = sp.Poly(expr, v, z)
        return sum(c * v**m[0] * z**m[1] for m, c in Pl.terms() if sum(m) == k)
    X3, Y3 = parte_grado(vdot, 3), parte_grado(zdot, 3)
    C3 = sp.expand(v * Y3 - z * X3)
    comprobar(sp.expand(C3 - z**4 / 2) == 0, f"C_3 = v Y_3 - z X_3 = {C3} = z^4/2")

    # primer blow-up z = v zeta, dividiendo por v^2
    Xv, Xz, k = explotar_campo(vdot, zdot, v, z, 1, v, zeta)
    print(f"\n    Tras z = v*zeta se divide por v^{k}:")
    print("    v'    =", Xv)
    print("    zeta' =", Xz)
    comprobar(k == 2, "la potencia común de v que se divide es v^2 (= v^(k-1), k = 3)")
    v_doc = (v * (P * R * zeta + P * Q * R * zeta**2 - zeta**3 / 2)
             + v**2 * (-sp.Rational(1, 2) - Q * zeta / 2 + (P * (S - 1) + sp.Rational(1, 2)) * zeta**2
                       + (P * Q * S + Q / 2) * zeta**3))
    z_doc = zeta**4 / 2 + v * (zeta / 2 + Q * zeta**2 / 2 - zeta**3 / 2 - Q * zeta**4 / 2)
    comprobar(sp.expand(Xv - v_doc) == 0, "v' coincide con la monografía")
    comprobar(sp.expand(Xz - z_doc) == 0, "zeta' coincide con la monografía")
    J1 = jac(Xv, Xz, v, zeta).subs({v: 0, zeta: 0})
    comprobar(J1 == sp.zeros(2, 2), "(v,zeta) = (0,0): parte lineal nula")
    def parte_grado2(expr, k):
        Pl = sp.Poly(expr, v, zeta)
        return sum(c * v**m[0] * zeta**m[1] for m, c in Pl.terms() if sum(m) == k)
    X2, Y2 = parte_grado2(Xv, 2), parte_grado2(Xz, 2)
    comprobar(sp.expand(X2 - (-v**2 / 2 + P * R * v * zeta)) == 0 and sp.expand(Y2 - v * zeta / 2) == 0,
              "parte cuadrática (-v^2/2 + PR v zeta, v zeta/2)")
    C2 = sp.factor(sp.expand(v * Y2 - zeta * X2))
    comprobar(sp.expand(C2 - v * zeta * (v - P * R * zeta)) == 0, f"C_2 = {C2} = v zeta (v - PR zeta)")

    # segundo blow-up en cada dirección: valores propios
    print("\n    Segundos blow-ups (tabla de valores propios):")
    def equilibrios_en_E(Xu, Xs, var_u, var_s):
        """Raíces de Xs(u=0, s) = 0 -> puntos de E, y valores propios allí."""
        pts = sp.solve(sp.Eq(Xs.subs(var_u, 0), 0), var_s)
        res = []
        for p in pts:
            Jp = jac(Xu, Xs, var_u, var_s).subs({var_u: 0, var_s: p})
            res.append((p, sp.simplify(Jp), {sp.simplify(k): m for k, m in Jp.eigenvals().items()}))
        return res
    # (A) carta zeta = v eta  (direccion zeta = 0)
    Xv2, Xe2, k2 = explotar_campo(Xv, Xz, v, zeta, 1, v, eta)
    print(f"    carta zeta = v eta, se divide por v^{k2}")
    for p, Jp, ev in equilibrios_en_E(Xv2, Xe2, v, eta):
        print(f"      punto eta = {p}: valores propios {ev}")
    # (B) carta v = zeta xi (direcciones v = 0 y v = PR zeta)
    Xz3, Xx3, k3 = explotar_campo(Xz, Xv, zeta, v, 1, zeta, xi)
    print(f"    carta v = zeta xi, se divide por zeta^{k3}")
    for p, Jp, ev in equilibrios_en_E(Xz3, Xx3, zeta, xi):
        print(f"      punto xi = {p}: valores propios {ev}")
    # comprobaciones explícitas de la tabla
    ev_A = equilibrios_en_E(Xv2, Xe2, v, eta)
    comprobar(any(p == 0 and set(ev) == {-sp.Rational(1, 2), 1} for p, _, ev in ev_A),
              "dirección zeta = 0: valores propios {-1/2, 1}  (silla)")
    ev_B = equilibrios_en_E(Xz3, Xx3, zeta, xi)
    comprobar(any(sp.simplify(p - P * R) == 0 and set(ev) == {-P * R, P * R / 2} for p, _, ev in ev_B),
              "dirección v = PR zeta: valores propios {-PR, PR/2}  (silla)")
    comprobar(any(p == 0 and set(ev) == {P * R, 0} for p, _, ev in ev_B),
              "dirección v = 0: valores propios {PR, 0}  (semihiperbólico)")
    # sobre xi = 0 : zeta' = zeta^3/2  (tras dividir)
    sobre = sp.expand(Xz3.subs(xi, 0))
    print(f"    sobre xi = 0:  zeta' = {sobre} + ...")
    comprobar(sp.expand(sobre - zeta**3 / 2) == 0, "zeta' = zeta^3/2 sobre xi = 0")
    # ---------- (d) pliegues para Q = 2 ----------
    print("\n(d) Pliegues de f para Q = 2 (régimen lento-rápido):")
    f2 = f.subs(Q, 2)
    ceros = sp.solve(sp.diff(f2, x), x)
    esperados = {-sp.Rational(1, 2), (-1 + sp.sqrt(5)) / 2, (-1 - sp.sqrt(5)) / 2}
    comprobar({sp.simplify(c) for c in ceros} == esperados, f"ceros de f' (Q=2): {ceros} = {{-1/2, (-1±√5)/2}}")
    comprobar(sp.discriminant(sp.diff(f2, x), x) != 0, "los tres ceros de f' son simples (discriminante no nulo)")
    comprobar([sp.simplify(f2.subs(x, c)) for c in [-sp.Rational(1, 2), (-1 + sp.sqrt(5)) / 2, (-1 - sp.sqrt(5)) / 2]]
              == [-sp.Rational(9, 16), 1, 1], "f en los pliegues: -9/16, 1, 1")
    Phi2 = (x / (x + Q) - R * f).subs(Q, 2)
    Sc = [sp.simplify(Phi2.subs(x, c)) for c in [-sp.Rational(1, 2), (-1 + sp.sqrt(5)) / 2, (-1 - sp.sqrt(5)) / 2]]
    comprobar(sp.simplify(Sc[0] - (9 * R / 16 - sp.Rational(1, 3))) == 0, "S_c(-1/2) = 9R/16 - 1/3")
    comprobar(sp.simplify(Sc[1] - (sp.sqrt(5) - 2 - R)) == 0, "S_c((-1+√5)/2) = √5 - 2 - R")
    comprobar(sp.simplify(Sc[2] - (-2 - sp.sqrt(5) - R)) == 0, "S_c((-1-√5)/2) = -2 - √5 - R < 0")
    print("\nOK: verificación simbólica del sistema completada.")
