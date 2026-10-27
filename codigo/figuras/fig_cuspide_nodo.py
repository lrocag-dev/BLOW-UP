"""Figura: cúspide y nodo con su transformada estricta en la carta (x,t), y = x t."""
import numpy as np

from estilo import plt, guardar, C_CURVA, C_E1

if __name__ == "__main__":
    fig, axs = plt.subplots(2, 2, figsize=(8.2, 7.0), constrained_layout=True)
    s = np.linspace(-1.3, 1.3, 600)

    # --- cúspide ---
    ax = axs[0, 0]
    ax.plot(s**2, s**3, color=C_CURVA, lw=1.6, label="$C:\\ y^2=x^3$")
    ax.plot(0, 0, "ko", ms=4)
    ax.set_title("Cúspide en el plano $(x,y)$")
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_aspect("equal"); ax.legend(loc="lower right")
    ax = axs[0, 1]
    t = np.linspace(-1.3, 1.3, 400)
    ax.plot(t**2, t, color=C_CURVA, lw=1.6, label="$\\tilde C:\\ t^2=x$")
    ax.axvline(0, color=C_E1, lw=2.0, label="$E:\\ x=0$")
    ax.plot(0, 0, "ko", ms=4)
    ax.annotate("$\\tilde C$ es lisa y\ntangente a $E$", (0, 0), xytext=(0.75, -0.9), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.set_title("Transformada estricta en la carta $(x,t)$, $y=xt$")
    ax.set_xlabel("$x$"); ax.set_ylabel("$t$"); ax.set_xlim(-0.5, 1.7); ax.set_ylim(-1.4, 1.4)
    ax.legend(loc="upper right")

    # --- nodo ---
    ax = axs[1, 0]
    x = np.linspace(-1.0, 1.2, 600)
    yy = x * np.sqrt(1 + x)
    ax.plot(x, yy, color=C_CURVA, lw=1.6, label="$C:\\ y^2=x^2(x+1)$")
    ax.plot(x, -yy, color=C_CURVA, lw=1.6)
    ax.plot(0, 0, "ko", ms=4)
    ax.set_title("Nodo en el plano $(x,y)$")
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_aspect("equal"); ax.legend(loc="lower right")
    ax = axs[1, 1]
    t = np.linspace(-1.6, 1.6, 400)
    ax.plot(t**2 - 1, t, color=C_CURVA, lw=1.6, label="$\\tilde C:\\ t^2=x+1$")
    ax.axvline(0, color=C_E1, lw=2.0, label="$E:\\ x=0$")
    ax.plot([0, 0], [-1, 1], "ko", ms=4)
    ax.annotate("$t=1$", (0, 1), xytext=(0.25, 1.2), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.annotate("$t=-1$", (0, -1), xytext=(0.25, -1.4), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.set_title("Transformada estricta en la carta $(x,t)$: dos puntos sobre $E$")
    ax.set_xlabel("$x$"); ax.set_ylabel("$t$"); ax.set_xlim(-1.2, 1.5); ax.set_ylim(-1.7, 1.7)
    ax.legend(loc="upper right")
    guardar(fig, "cuspide_nodo.pdf")
