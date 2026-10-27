"""Blow-up direccional del origen del sistema de la monografía y simulación numérica.

  carta 1:  x = r,  y = r v      carta 2:  y = r,  x = r u        (E = {r = 0})

Para k = 1 (parte lineal no nula) el campo en la carta ya es regular, sin dividir.
Se integra con scipy.integrate.solve_ivp (DOP853, rtol=1e-10, atol=1e-12, semilla fija)
y se genera el retrato de fase antes y después del blow-up.

EVIDENCIA NUMÉRICA, NO DEMOSTRACIÓN.  Parámetros: conjunto A (P=1, Q=0.3, R=1, S=0.33).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from sistema import (campo, campo_carta1, campo_carta2, f, g, jacobiana, equilibrios, clasificar,
                     RTOL, ATOL, METODO, SEMILLA)

SALIDA_FIG = Path(__file__).resolve().parents[1] / "figuras" / "figuras_salida"
SALIDA_FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"pdf.fonttype": 42, "font.size": 9, "axes.grid": True, "grid.alpha": 0.25})

P, Q, R, S = 1.0, 0.3, 1.0, 0.33
ARGS = (P, Q, R, S)
V_ESTRELLA = Q + 2 * P * Q * S          # equilibrio de E en la carta 1 (v = y/x de la dirección estable)
U_ESTRELLA = 1.0 / V_ESTRELLA           # equilibrio de E en la carta 2


def integrar(rhs, y0, T, ventana, n=2000):
    """Integra de 0 a T (T<0: hacia atrás) hasta salir de ventana=(xmin,xmax,ymin,ymax)."""
    xmin, xmax, ymin, ymax = ventana

    def salida(t, s, *a):
        return min(xmax - s[0], s[0] - xmin, ymax - s[1], s[1] - ymin)
    salida.terminal, salida.direction = True, -1
    sol = solve_ivp(rhs, (0.0, T), y0, args=ARGS, method=METODO, rtol=RTOL, atol=ATOL,
                    events=salida, t_eval=np.linspace(0.0, T, n))
    return sol


def trayectoria(rhs, y0, T, ventana):
    f1 = integrar(rhs, y0, T, ventana)
    b1 = integrar(rhs, y0, -T, ventana)
    return np.hstack([b1.y[:, ::-1], f1.y])


def flecha(ax, XY, color):
    """Flecha en el punto medio de la órbita medido por longitud de arco."""
    if XY.shape[1] < 10:
        return
    d = np.hypot(np.diff(XY[0]), np.diff(XY[1]))
    acum = np.concatenate([[0.0], np.cumsum(d)])
    i = int(np.searchsorted(acum, acum[-1] / 2))
    i = min(max(i, 0), XY.shape[1] - 2)
    ax.annotate("", xy=(XY[0, i + 1], XY[1, i + 1]), xytext=(XY[0, i], XY[1, i]),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8, mutation_scale=8))


def dibujar_trayectorias(ax, rhs, ics, T, ventana, color="tab:blue", pre=None):
    for y0 in ics:
        XY = trayectoria(rhs, y0, T, ventana)
        if pre is not None:
            XY = pre(XY)
        ax.plot(XY[0], XY[1], color=color, lw=0.7)
        flecha(ax, XY, color)


def marcar_equilibrios(ax, eq, filtro=None):
    for (x0, y0, _) in eq:
        if filtro is not None and not filtro(x0, y0):
            continue
        J = jacobiana(x0, y0, *ARGS)
        tipo, hip, lam = clasificar(J)
        if "silla" in tipo:
            ax.plot(x0, y0, marker="D", mfc="white", mec="k", ms=6, ls="none")
        elif "inestable" in tipo:
            ax.plot(x0, y0, marker="o", mfc="white", mec="k", ms=6, ls="none")
        elif "estable" in tipo:
            ax.plot(x0, y0, marker="o", mfc="k", mec="k", ms=6, ls="none")
        else:
            ax.plot(x0, y0, marker="s", mfc="tab:red", mec="k", ms=6, ls="none")


def leyenda_equilibrios(ax):
    ax.plot([], [], marker="D", mfc="white", mec="k", ls="none", label="silla")
    ax.plot([], [], marker="o", mfc="white", mec="k", ls="none", label="inestable")
    ax.plot([], [], marker="o", mfc="k", mec="k", ls="none", label="estable")


if __name__ == "__main__":
    rng = np.random.default_rng(SEMILLA)
    print("=" * 78)
    print("BLOW-UP DIRECCIONAL DEL ORIGEN: simulación numérica (evidencia, no demostración)")
    print(f"solve_ivp: método={METODO}, rtol={RTOL}, atol={ATOL}; semilla={SEMILLA}")
    print(f"parámetros: P={P}, Q={Q}, R={R}, S={S}")
    print("=" * 78)

    # ---- 1. Equilibrios de E en las cartas y su linealización ----
    print(f"\nEquilibrios sobre E=\\{{r=0\\}}:  carta 1: v = 0 y v* = Q+2PQS = {V_ESTRELLA:.6f};  carta 2: u* = 1/v* = {U_ESTRELLA:.6f}")
    for nombre, rhs, pts in [("carta 1", campo_carta1, [(0.0, 0.0), (0.0, V_ESTRELLA)]),
                             ("carta 2", campo_carta2, [(0.0, U_ESTRELLA)])]:
        for p in pts:
            res = np.max(np.abs(rhs(0.0, p, *ARGS)))
            assert res < 1e-12, f"{nombre}: {p} no es equilibrio (residuo {res})"
            h = 1e-6
            Jn = np.column_stack([(np.array(rhs(0, [p[0] + h * (k == 0), p[1] + h * (k == 1)], *ARGS))
                                   - np.array(rhs(0, [p[0] - h * (k == 0), p[1] - h * (k == 1)], *ARGS))) / (2 * h)
                                  for k in range(2)])
            lam = np.sort(np.linalg.eigvals(Jn).real)
            print(f"  {nombre}: (r,{'v' if nombre=='carta 1' else 'u'}) = ({p[0]:.4f}, {p[1]:.4f})  valores propios ≈ {np.round(lam, 6)}")
    # predicción analítica: carta 1 en (0,0): Q/2, -(PQS+Q/2);  en (0,v*): -PQS, Q/2+PQS
    #                      (los valores propios de la jacobiana del campo original en el origen son Q/2 y -PQS)
    esperados1 = {(0.0, 0.0): sorted([Q / 2, -(P * Q * S + Q / 2)]), (0.0, V_ESTRELLA): sorted([-P * Q * S, Q / 2 + P * Q * S])}
    for p, esp in esperados1.items():
        h = 1e-6
        Jn = np.column_stack([(np.array(campo_carta1(0, [p[0] + h * (k == 0), p[1] + h * (k == 1)], *ARGS))
                               - np.array(campo_carta1(0, [p[0] - h * (k == 0), p[1] - h * (k == 1)], *ARGS))) / (2 * h)
                              for k in range(2)])
        lam = np.sort(np.linalg.eigvals(Jn).real)
        assert np.allclose(lam, esp, atol=1e-7), f"valores propios en carta 1 {p}: {lam} vs {esp}"
    print("  coinciden con las predicciones analíticas (carta 1): sillas hiperbólicas sobre E.")

    # ---- 2. Validación cruzada: órbitas de las cartas vs órbita directa ----
    print("\nValidación cruzada (misma parametrización del tiempo, k=1 no reescala):")
    T_COMP = 20.0
    XMIN, YMIN = 1e-3, 1e-5   # las cartas dejan de ser válidas si x -> 0 (carta 1) o y -> 0 (carta 2)

    def sale_de_carta(t, s, *a):
        return min(s[0] - XMIN, s[1] - YMIN)
    sale_de_carta.terminal, sale_de_carta.direction = True, -1

    desv1, desv2 = 0.0, 0.0
    n_ic = 8
    for k in range(n_ic):
        r0, v0 = rng.uniform(0.01, 0.2), rng.uniform(0.05, 2.5)
        x0, y0 = r0, r0 * v0
        directa = solve_ivp(campo, (0.0, T_COMP), [x0, y0], args=ARGS, method=METODO, rtol=RTOL, atol=ATOL,
                            events=sale_de_carta)
        T_fin = directa.t[-1]
        tt = np.linspace(0.0, T_fin, 401)
        directa = solve_ivp(campo, (0.0, T_fin), [x0, y0], args=ARGS, method=METODO, rtol=RTOL, atol=ATOL, t_eval=tt)
        c1 = solve_ivp(campo_carta1, (0.0, T_fin), [r0, v0], args=ARGS, method=METODO, rtol=RTOL, atol=ATOL, t_eval=tt)
        err = np.max(np.abs(np.vstack([c1.y[0], c1.y[0] * c1.y[1]]) - directa.y))
        desv1 = max(desv1, err)
        u0 = x0 / y0   # carta 2 (requiere y0 > 0)
        c2 = solve_ivp(campo_carta2, (0.0, T_fin), [y0, u0], args=ARGS, method=METODO, rtol=RTOL, atol=ATOL, t_eval=tt)
        err2 = np.max(np.abs(np.vstack([c2.y[0] * c2.y[1], c2.y[0]]) - directa.y))
        desv2 = max(desv2, err2)
        print(f"  IC {k+1}: (x0,y0)=({x0:.4f},{y0:.4f}), T_comparación={T_fin:6.3f}  "
              f"max|π(carta 1)-directa|={err:.2e}  max|π(carta 2)-directa|={err2:.2e}")
        assert directa.success and c1.success and c2.success
    print(f"  desviación máxima: carta 1 = {desv1:.2e}, carta 2 = {desv2:.2e}")
    assert desv1 < 1e-6 and desv2 < 1e-6, "las órbitas en las cartas no coinciden con la órbita directa"

    # ---- 3. Retrato de fase ANTES ----
    eq, _ = equilibrios(*ARGS)
    eq1 = [e for e in eq if e[0] >= -1e-12 and e[1] >= -1e-12]
    rng = np.random.default_rng(SEMILLA)
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    # (a) global
    ax = axs[0]
    vent = (-0.02, 1.15, -0.02, 0.5)
    ics = [(rng.uniform(0.02, 1.1), rng.uniform(0.01, 0.45)) for _ in range(40)] + [(0.5, 0.0), (0.2, 0.0), (0.9, 0.0)]
    dibujar_trayectorias(ax, campo, ics, 60.0, (-0.02, 1.15, -0.0, 0.5))
    xs = np.linspace(0, 1.15, 400)
    ax.plot(xs, f(xs, Q), "k--", lw=0.9, label="$\\dot x=0$: $y=f(x)$")
    ax.plot(xs, (xs / (xs + Q) - S) / R, "k:", lw=0.9, label="$g=0$")
    marcar_equilibrios(ax, eq1)
    leyenda_equilibrios(ax)
    ax.set_xlim(vent[0], vent[1]); ax.set_ylim(-0.02, 0.5)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_title("(a) Antes del blow-up: primer cuadrante")
    ax.legend(loc="upper right", fontsize=7)
    # (b) zoom
    ax = axs[1]
    vz = (-0.005, 0.12, -0.003, 0.06)
    ics_z = [(rng.uniform(0.002, 0.1), rng.uniform(0.0005, 0.05)) for _ in range(40)]
    dibujar_trayectorias(ax, campo, ics_z, 120.0, vz)
    ax.set_xlim(vz[0], vz[1]); ax.set_ylim(vz[2], vz[3])
    marcar_equilibrios(ax, eq1, lambda a, b: a < 0.12 and b < 0.06)
    leyenda_equilibrios(ax); ax.legend(loc="lower right", fontsize=7)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_title("(b) Entorno del origen (silla)")
    fig.savefig(SALIDA_FIG / "dinamica_antes.pdf", metadata={"CreationDate": None})
    plt.close(fig)

    # ---- 4. Retrato de fase DESPUÉS ----
    rng = np.random.default_rng(SEMILLA + 1)
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.2), constrained_layout=True)
    # carta 1
    ax = axs[0]
    v1 = (-0.01, 0.30, -0.05, 2.5)
    ics1 = [(rng.uniform(0.005, 0.28), rng.uniform(0.02, 2.4)) for _ in range(45)]
    ics1 += [(0.1, 0.0), (0.2, 0.0), (0.1, V_ESTRELLA), (0.2, V_ESTRELLA)]
    dibujar_trayectorias(ax, campo_carta1, ics1, 120.0, (0.0, v1[1], v1[2] + 0.04, v1[3]), color="tab:green")
    ax.axvline(0, color="k", lw=2.0)
    ax.plot([0, 0], [0, V_ESTRELLA], ls="none", marker="D", mfc="white", mec="k", ms=6)
    ax.set_xlim(v1[0], v1[1]); ax.set_ylim(v1[2], v1[3])
    ax.set_xlabel("$r$"); ax.set_ylabel("$v=y/x$")
    ax.set_title("(a) Carta 1: $x=r,\\ y=rv$")
    ax.text(0.012, 2.3, "$E$ ($r=0$)", fontsize=8, bbox=dict(fc="white", ec="none", alpha=0.85))
    # carta 2
    ax = axs[1]
    v2 = (-0.01, 0.30, -0.3, 8.0)
    ics2 = [(rng.uniform(0.005, 0.28), rng.uniform(0.05, 7.5)) for _ in range(45)]
    ics2 += [(0.1, U_ESTRELLA), (0.2, U_ESTRELLA)]
    dibujar_trayectorias(ax, campo_carta2, ics2, 120.0, (0.0, v2[1], v2[2] + 0.2, v2[3]), color="tab:purple")
    ax.axvline(0, color="k", lw=2.0)
    ax.plot([0], [U_ESTRELLA], ls="none", marker="D", mfc="white", mec="k", ms=6)
    ax.set_xlim(v2[0], v2[1]); ax.set_ylim(v2[2], v2[3])
    ax.set_xlabel("$r$"); ax.set_ylabel("$u=x/y$")
    ax.set_title("(b) Carta 2: $y=r,\\ x=ru$")
    # imagen por pi de la carta 1 sobre las órbitas directas
    ax = axs[2]
    vz = (-0.005, 0.12, -0.003, 0.06)
    for y0 in [(r0, r0 * v0) for r0, v0 in [(0.03, 0.2), (0.05, 1.0), (0.08, 0.4), (0.02, 0.7), (0.1, 0.1)]]:
        XYd = trayectoria(campo, list(y0), 60.0, (-0.02, 0.2, -0.01, 0.1))
        ax.plot(XYd[0], XYd[1], color="tab:blue", lw=2.2, alpha=0.5)
        r0, v0 = y0[0], y0[1] / y0[0]
        XYc = trayectoria(campo_carta1, [r0, v0], 60.0, (0.0, 0.2, 0.0, 3.0))
        ax.plot(XYc[0], XYc[0] * XYc[1], color="tab:green", lw=0.8, ls="--")
    ax.plot([], [], color="tab:blue", lw=2.2, alpha=0.5, label="órbita directa $(x,y)$")
    ax.plot([], [], color="tab:green", lw=0.8, ls="--", label="$\\pi$(órbita en la carta 1)")
    ax.set_xlim(vz[0], vz[1]); ax.set_ylim(vz[2], vz[3])
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_title("(c) $\\pi$ de las órbitas de (a) vs. directas")
    ax.legend(fontsize=7, loc="upper left", framealpha=0.9)
    fig.savefig(SALIDA_FIG / "dinamica_despues.pdf", metadata={"CreationDate": None})
    plt.close(fig)
    print(f"\nFiguras escritas en {SALIDA_FIG}: dinamica_antes.pdf, dinamica_despues.pdf")
    print("OK: simulación en las cartas coherente con la dinámica original (evidencia numérica).")
