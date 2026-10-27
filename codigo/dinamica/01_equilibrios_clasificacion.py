"""Localización numérica de equilibrios en el primer cuadrante, jacobiano y
clasificación (marcando los no hiperbólicos), para varias elecciones de parámetros.

EVIDENCIA NUMÉRICA, NO DEMOSTRACIÓN: aritmética de punto flotante.

Método: raíces reales de h(x) (polinomio de grado 5, numpy.roots + Newton) y
contraste independiente con búsqueda multi-arranque (scipy.optimize.fsolve, semilla fija).
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve

from sistema import (f, fp, g, campo, jacobiana, equilibrios, clasificar, polinomio_h,
                     SEMILLA, TOL_HIPERBOLICO, TOL_RESIDUO)

SALIDA = Path(__file__).resolve().parent / "salida"
SALIDA.mkdir(exist_ok=True)

# (nombre, P, Q, R, S)   Q = 0.3, R = 1 como en la monografía salvo indicación
CONJUNTOS = [
    ("A", 1.0, 0.3, 1.0, 0.33),   # 3 equilibrios interiores
    ("B", 1.0, 0.3, 1.0, 0.36),   # 1 interior
    ("C", 1.0, 0.3, 1.0, 1.0),    # 1 interior
    ("D", 1.0, 0.3, 1.0, 1.5),    # raíces con x>1: y<0, fuera del cuadrante
    ("E", 1.0, 0.5, 1.0, 2.0 / 3.0),  # S(1+Q)=1: (1,0) semihiperbólico
]


def en_primer_cuadrante(x, y, tol=1e-12):
    return x >= -tol and y >= -tol


def multiarranque(P, Q, R, S, n=400):
    """Equilibrios del primer cuadrante hallados con fsolve desde arranques aleatorios (semilla fija)."""
    rng = np.random.default_rng(SEMILLA)
    hallados = []
    for _ in range(n):
        x0, y0 = rng.uniform(0.0, 1.3), rng.uniform(0.0, 0.8)
        sol, info, ier, _ = fsolve(lambda s: campo(0.0, s, P, Q, R, S), [x0, y0], full_output=True, xtol=1e-14)
        if ier == 1 and np.max(np.abs(campo(0.0, sol, P, Q, R, S))) < TOL_RESIDUO and en_primer_cuadrante(*sol, tol=1e-9):
            if not any(np.hypot(*(sol - w)) < 1e-7 for w in hallados):
                hallados.append(sol)
    return hallados


def resolver_bogdanov_takens(Q=0.3, R=1.0, x_ini=0.575):
    """Resuelve det J = tr J = 0 en (x*, P, S): Phi'(x*) = 0, f'(x*) = 2 P R (x*+Q) f(x*), S = Phi(x*)."""
    from numpy.polynomial import Polynomial as Pl
    pf = Pl([0.0, Q, 1.0, -Q, -1.0])
    # Phi'(x) = Q/(x+Q)^2 - R f'(x)  =>  Q - R f'(x) (x+Q)^2 = 0
    num = Pl([Q]) - R * pf.deriv() * Pl([Q, 1.0]) ** 2
    raices = [r.real for r in num.roots() if abs(r.imag) < 1e-9]
    xs = min(raices, key=lambda r: abs(r - x_ini))
    for _ in range(30):  # Newton
        xs -= num(xs) / num.deriv()(xs)
    S = xs / (xs + Q) - R * f(xs, Q)
    P = fp(xs, Q) / (2.0 * R * (xs + Q) * f(xs, Q))
    return xs, f(xs, Q), P, S


if __name__ == "__main__":
    filas_tex = []
    print("=" * 78)
    print("EQUILIBRIOS EN EL PRIMER CUADRANTE (evidencia numérica, punto flotante)")
    print("(* = raíz múltiple de h agrupada: el punto flotante no distingue una raíz doble de dos simples)")
    print(f"semilla = {SEMILLA};  tolerancia de no-hiperbolicidad = {TOL_HIPERBOLICO};  residuo máx. = {TOL_RESIDUO}")
    print("=" * 78)
    todos = list(CONJUNTOS)
    # Conjunto F: candidato Bogdanov-Takens de la monografía (Q = 0.3, R = 1)
    xs_bt, ys_bt, P_bt, S_bt = resolver_bogdanov_takens()
    todos.append(("F", P_bt, 0.3, 1.0, S_bt))
    for nombre, P, Q, R, S in todos:
        print(f"\n--- Conjunto {nombre}: P = {P:.6g}, Q = {Q}, R = {R}, S = {S:.6g} ---")
        eq, dudosas = equilibrios(P, Q, R, S)
        if dudosas:
            print(f"    AVISO: raíces de h casi reales (|Im| < 1e-4 pero > {1e-9}): {dudosas}")
        propios = [e for e in eq if en_primer_cuadrante(e[0], e[1])]
        for (x0, y0, lugar) in propios:
            res = np.max(np.abs(campo(0.0, [x0, y0], P, Q, R, S)))
            assert res < (1e-7 if lugar.endswith('*') else TOL_RESIDUO), f"residuo demasiado grande en ({x0},{y0}): {res}"
            J = jacobiana(x0, y0, P, Q, R, S)
            tipo, hip, lam = clasificar(J)
            marca = "" if hip else "   <-- NO HIPERBÓLICO"
            ls = ", ".join(f"{l.real:+.4f}{l.imag:+.4f}i" if abs(l.imag) > 1e-12 else f"{l.real:+.4f}" for l in lam)
            print(f"    ({x0:+.6f}, {y0:+.6f}) [{lugar:9s}]  det={np.linalg.det(J):+.3e} tr={np.trace(J):+.3e}  "
                  f"lambda = {ls}  ->  {tipo}{marca}")
            filas_tex.append((nombre, x0, y0, lam, tipo, hip, lugar))
        # contraste independiente
        mult = multiarranque(P, Q, R, S)
        for w in mult:
            assert any(np.hypot(w[0] - e[0], w[1] - e[1]) < 1e-6 for e in propios), \
                f"fsolve halló un equilibrio {w} que no está en la lista de raíces"
        for e in propios:
            if e[1] > 1e-9 and e[0] > 0:  # interiores
                assert any(np.hypot(w[0] - e[0], w[1] - e[1]) < 1e-6 for w in mult), \
                    f"no se halló por multi-arranque el equilibrio interior {e}"
        print(f"    contraste fsolve multi-arranque: {len(mult)} equilibrios del primer cuadrante, coinciden con las raíces de h")

    # ---------- afirmaciones numéricas de la monografía ----------
    print("\n" + "=" * 78)
    print("COMPROBACIÓN NUMÉRICA DE AFIRMACIONES NUMÉRICAS DE LA MONOGRAFÍA (Q = 0.3, R = 1)")
    print("=" * 78)
    Q, R = 0.3, 1.0
    from numpy.polynomial import Polynomial as Pl
    pf = Pl([0.0, Q, 1.0, -Q, -1.0])
    num = Pl([Q]) - R * pf.deriv() * Pl([Q, 1.0]) ** 2
    crit = sorted(float(c.real) for c in num.roots() if abs(c.imag) < 1e-9)
    Phi = lambda x: x / (x + Q) - R * f(x, Q)
    valores = sorted(Phi(c) for c in crit)
    print(f"  valores críticos de Phi (raíces reales de Phi'): {[round(float(v), 4) for v in valores]}")
    print(f"  monografía afirma: ≈ 0.3204, 0.3388, 1.3962")
    esperado = [0.3204, 0.3388, 1.3962]
    # Los valores críticos con x > -Q.  (Phi' tiene además raíces con x < -Q)
    mios = [round(float(v), 4) for v in valores]
    print(f"  (todas las raíces reales de Phi' dan los valores {mios})")
    for e in esperado:
        assert any(abs(v - e) < 6e-4 for v in valores), f"valor crítico {e} no reproducido"
    for S, cuantas in [(0.33, 3), (0.36, 1), (1.0, 1), (1.5, 3)]:
        raices, _ = __import__("sistema").raices_reales(polinomio_h(1.0, Q, R, S))
        print(f"  S = {S}: raíces reales de h = {[round(float(r), 3) for r in raices]}  ({len(raices)})")
        assert len(raices) == cuantas, f"S={S}: se esperaban {cuantas} raíces reales"
    raices, _ = __import__("sistema").raices_reales(polinomio_h(1.0, Q, R, 0.33))
    for r, e in zip(raices, [0.266, 0.457, 0.660]):
        assert abs(r - e) < 6e-4, f"raíz {r} difiere de la citada {e}"
    print(f"  Bogdanov-Takens: x* = {xs_bt:.4f}, y* = {ys_bt:.4f}, P = {P_bt:.4f}, S = {S_bt:.4f} "
          f"(monografía: 0.5751, 0.3368, 0.6646, 0.3204)")
    for val, ref in [(xs_bt, 0.5751), (ys_bt, 0.3368), (P_bt, 0.6646), (S_bt, 0.3204)]:
        assert abs(val - ref) < 6e-4, f"valor BT {val} difiere de {ref}"
    J = jacobiana(xs_bt, ys_bt, P_bt, 0.3, 1.0, S_bt)
    print(f"  J = {np.round(J, 4).tolist()},  det J = {np.linalg.det(J):.2e},  tr J = {np.trace(J):.2e}")
    assert abs(np.linalg.det(J)) < 1e-12 and abs(np.trace(J)) < 1e-12 and np.abs(J).max() > 0.1
    # Fragmento LaTeX con la tabla (lo incluye el apéndice con \\input, así los valores coinciden)
    params = {n: (P, Q, R, S) for n, P, Q, R, S in todos}
    # El fragmento NO termina en \\ (lo añade el documento antes de \\bottomrule).
    lineas = []
    previo = None
    for (nombre, x0, y0, lam, tipo, hip, lugar) in filas_tex:
        if nombre != previo:
            if previo is not None:
                lineas.append("\\midrule")
            P_, Q_, R_, S_ = params[nombre]
            lineas.append(f"\\emph{{Conjunto {nombre}}} & \\multicolumn{{3}}{{l}}{{$P={P_:.4f},\\ Q={Q_:g},\\ R={R_:g},\\ S={S_:.4f}$}} \\\\")
            previo = nombre
        if abs(lam[0].imag) > 1e-12:
            ls = f"${lam[0].real:+.4f}\\pm{abs(lam[0].imag):.4f}\\,i$"
        else:
            ls = "$" + ",\\ ".join(f"{l.real:+.4f}" for l in lam) + "$"
        marca = "" if hip else "\\ \\ensuremath{{}^{\\dagger}}"
        lineas.append(f"\\quad & $({x0:.4f},\\ {y0:.4f})$ & {ls} & {tipo.split(' (')[0]}{marca} \\\\")
    lineas[-1] = lineas[-1][:-3].rstrip()
    with open(SALIDA / "tabla_equilibrios.tex", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lineas) + "\n")
    print(f"Fragmento LaTeX escrito en {SALIDA / 'tabla_equilibrios.tex'}")
    print("OK: equilibrios, jacobianos y clasificación calculados; contrastes superados.")
