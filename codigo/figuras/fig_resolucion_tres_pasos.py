"""Figura: resolución embebida de y^2 = x^3 en tres blow-ups.

Cartas usadas (como en el script simbolico/05_resolucion_cuspide.py):
  paso 0: (x,y)         curva C
  paso 1: (x,t),  y = x t     C~: t^2 = x (tangente a E1: x = 0)
  paso 2: (t,w),  x = t w     C~: t = w,  E1: w = 0,  E2: t = 0   (tres rectas concurrentes)
  paso 3: (t,s),  w = t s     C~: s = 1,  E1: s = 0,  E3: t = 0   (SNC; E2 está en la otra carta, s = infinito)
"""
import numpy as np

from estilo import plt, guardar, C_CURVA, C_E1, C_E2, C_E3

if __name__ == "__main__":
    fig, axs = plt.subplots(1, 4, figsize=(13.5, 3.6), constrained_layout=True)
    L = np.linspace(-1.5, 1.5, 400)

    ax = axs[0]
    s = np.linspace(-1.3, 1.3, 400)
    ax.plot(s**2, s**3, color=C_CURVA, lw=1.6, label="$C$")
    ax.plot(0, 0, "ko", ms=4)
    ax.set_xlim(-0.5, 1.7); ax.set_ylim(-1.6, 1.6)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_title("Paso 0: $y^2=x^3$\n(mult. 2 en el origen)")
    ax.legend(loc="lower right")

    ax = axs[1]
    t = np.linspace(-1.3, 1.3, 400)
    ax.plot(t**2, t, color=C_CURVA, lw=1.6, label="$\\tilde C$")
    ax.axvline(0, color=C_E1, lw=2.0, label="$E_1$ (coef. 2)")
    ax.plot(0, 0, "ko", ms=4)
    ax.set_xlim(-0.5, 1.7); ax.set_ylim(-1.6, 1.6)
    ax.set_xlabel("$x$"); ax.set_ylabel("$t$")
    ax.set_title("Paso 1: carta $y=xt$\n$\\tilde C$ y $E_1$ tangentes")
    ax.legend(loc="lower right")

    ax = axs[2]
    ax.plot(L, L, color=C_CURVA, lw=1.6, label="$\\tilde C:\\ t=w$")
    ax.axhline(0, color=C_E1, lw=2.0, label="$E_1:\\ w=0$")
    ax.axvline(0, color=C_E2, lw=2.0, label="$E_2:\\ t=0$ (coef. 3)")
    ax.plot(0, 0, "ko", ms=4)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("$t$"); ax.set_ylabel("$w$")
    ax.set_title("Paso 2: carta $x=tw$\ntres ramas concurrentes")
    ax.legend(loc="lower right", fontsize=7)

    ax = axs[3]
    ax.axhline(1, color=C_CURVA, lw=1.6, label="$\\tilde C:\\ s=1$")
    ax.axhline(0, color=C_E1, lw=2.0, label="$E_1:\\ s=0$")
    ax.axvline(0, color=C_E3, lw=2.0, label="$E_3:\\ t=0$ (coef. 6)")
    ax.plot([0, 0], [0, 1], "ko", ms=4)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-0.7, 1.7)
    ax.set_xlabel("$t$"); ax.set_ylabel("$s$")
    ax.set_title("Paso 3: carta $w=ts$\n$E_2$ está en la otra carta ($s=\\infty$)")
    ax.legend(loc="center right", fontsize=7)
    fig.suptitle("Resolución embebida de $y^2=x^3$:  $\\pi^*C=\\tilde C+2E_1+3E_2+6E_3$, divisor final de cruces normales",
                 fontsize=10)
    guardar(fig, "resolucion_cuspide_tres_pasos.pdf")
