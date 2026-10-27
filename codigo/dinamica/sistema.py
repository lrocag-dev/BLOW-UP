"""Sistema de la monografía (Cap. de dinámica, Sec. 'Aplicación al sistema propuesto'):

    x' = (f(x) - y)/2,   y' = P y g(x,y),
    f(x) = (1 - x^2) x (x + Q),   g(x,y) = x - (S + R y)(x + Q).

Aquí todo es NUMÉRICO (punto flotante).  Las identidades algebraicas exactas
están en 00_verificacion_simbolica.py.  Se incluyen también los campos en las
cartas del blow-up direccional del origen:

  carta 1:  x = r,   y = r v     (E = {r = 0})
  carta 2:  y = r,   x = r u     (E = {r = 0})
"""
import numpy as np
from numpy.polynomial import Polynomial

# ---- tolerancias y semillas (documentadas en el apéndice) ----
RTOL = 1e-10          # tolerancia relativa de solve_ivp
ATOL = 1e-12          # tolerancia absoluta de solve_ivp
METODO = "DOP853"     # Runge-Kutta explícito de orden 8
TOL_HIPERBOLICO = 1e-4   # |Re(lambda)| <= TOL => se marca como NO hiperbólico
TOL_RAIZ_MULTIPLE = 1e-5  # raíces de h a distancia < TOL se agrupan (raíz múltiple, mal condicionada)
TOL_RESIDUO = 1e-12      # residuo máximo admitido en un equilibrio refinado
TOL_IMAG = 1e-9          # |Im| <= TOL => raíz real
SEMILLA = 20261027       # semilla de numpy.random.default_rng


def f(x, Q):
    return (1.0 - x**2) * x * (x + Q)


def fp(x, Q):
    return (1.0 - 3.0 * x**2) * (x + Q) + (x - x**3)


def g(x, y, Q, R, S):
    return x - (S + R * y) * (x + Q)


def campo(t, s, P, Q, R, S):
    x, y = s
    return [0.5 * (f(x, Q) - y), P * y * g(x, y, Q, R, S)]


def jacobiana(x, y, P, Q, R, S):
    gg = g(x, y, Q, R, S)
    gx = 1.0 - (S + R * y)
    gy = -R * (x + Q)
    return np.array([[0.5 * fp(x, Q), -0.5],
                     [P * y * gx, P * gg + P * y * gy]])


# ---- cartas del blow-up del origen ----
def campo_carta1(t, s, P, Q, R, S):
    """x = r, y = r v:   r' = r((1-r^2)(r+Q) - v)/2,  v' = v[P g(r,rv) - ((1-r^2)(r+Q) - v)/2]."""
    r, v = s
    a = (1.0 - r**2) * (r + Q) - v
    return [0.5 * r * a, v * (P * g(r, r * v, Q, R, S) - 0.5 * a)]


def campo_carta2(t, s, P, Q, R, S):
    """y = r, x = r u:   r' = P r g(ru, r),  u' = (f(ru)/r - 1)/2 - u P g(ru, r)."""
    r, u = s
    gg = g(r * u, r, Q, R, S)
    return [P * r * gg, 0.5 * (u * (1.0 - (r * u) ** 2) * (r * u + Q) - 1.0) - u * P * gg]


# ---- equilibrios ----
def polinomio_h(P, Q, R, S):
    pf = Polynomial([0.0, Q, 1.0, -Q, -1.0])  # f(x) = -x^4 - Q x^3 + x^2 + Q x
    return R * Polynomial([Q, 1.0]) * pf + Polynomial([S * Q, -(1.0 - S)])


def raices_reales(poli, tol_imag=TOL_IMAG):
    """Raíces reales (refinadas con Newton) y lista de raíces dudosas (|Im| pequeño pero no nulo)."""
    rs = poli.roots()
    reales, dudosas = [], []
    dp = poli.deriv()
    for r in rs:
        if abs(r.imag) <= tol_imag:
            xr = r.real
            for _ in range(20):
                d = dp(xr)
                if d == 0:
                    break
                xr -= poli(xr) / d
            reales.append(xr)
        elif abs(r.imag) < 1e-4:
            dudosas.append(r)
    return sorted(reales), dudosas


def equilibrios(P, Q, R, S):
    """Todos los equilibrios reales: lista de (x, y, etiqueta)."""
    eq = [(x0, 0.0, "eje") for x0 in (0.0, 1.0, -1.0, -Q)]
    # Q = 1 => -1 y -Q coinciden: se elimina el duplicado
    vistos, res = [], []
    for e in eq:
        if not any(abs(e[0] - w) < 1e-12 for w in vistos):
            vistos.append(e[0])
            res.append(e)
    reales, dudosas = raices_reales(polinomio_h(P, Q, R, S))
    grupos = []
    for xs in reales:  # agrupar raíces casi coincidentes (raíz múltiple: el punto flotante no distingue)
        if grupos and abs(xs - grupos[-1][-1]) < TOL_RAIZ_MULTIPLE:
            grupos[-1].append(xs)
        else:
            grupos.append([xs])
    for grupo in grupos:
        xs = float(np.mean(grupo))
        ys = f(xs, Q)
        if abs(ys) > 1e-9:  # y = 0 ya está en el eje
            res.append((xs, ys, "interior*" if len(grupo) > 1 else "interior"))
    return res, dudosas


def clasificar(J, tol=TOL_HIPERBOLICO):
    """Devuelve (tipo, es_hiperbolico, valores propios)."""
    lam = np.linalg.eigvals(J)
    re, im = lam.real, lam.imag
    escala = max(1.0, np.abs(J).max())
    if np.abs(J).max() <= tol:
        return "totalmente degenerado (J = 0)", False, lam
    if np.all(np.abs(re) <= tol * escala):
        if np.all(np.abs(lam) <= tol * escala):
            return "nilpotente (ambos valores propios 0, J != 0)", False, lam
        return "centro lineal (parte real nula)", False, lam
    if np.any(np.abs(re) <= tol * escala):
        return "semihiperbólico (un valor propio 0)", False, lam
    if np.all(np.abs(im) <= tol * escala):
        if re[0] * re[1] < 0:
            return "silla", True, lam
        return ("nodo estable" if re[0] < 0 else "nodo inestable"), True, lam
    return ("foco estable" if re[0] < 0 else "foco inestable"), True, lam
