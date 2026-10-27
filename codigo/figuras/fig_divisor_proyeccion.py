"""Figura: divisor excepcional E y la proyección pi: (x,t) -> (x, x t)."""
import numpy as np

from estilo import plt, guardar, C_CURVA, C_E1

if __name__ == "__main__":
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.3), constrained_layout=True)
    cs = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]
    cmap = plt.get_cmap("viridis")
    cols = [cmap(i / (len(cs) - 1)) for i in range(len(cs))]
    xs = np.linspace(-1.5, 1.5, 200)

    ax = axs[0]
    for c, col in zip(cs, cols):
        ax.plot(xs, np.full_like(xs, c), color=col, lw=1.1, label=f"$t={c:g}$")
    ax.axvline(0, color=C_E1, lw=3.0, label="$E=\\{x=0\\}$")
    t = np.linspace(-1.4, 1.4, 300)
    ax.plot(t**2, t, color=C_CURVA, lw=1.6, label="$\\tilde C:\\ t^2=x$")
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-2.5, 2.5)
    ax.set_xlabel("$x$"); ax.set_ylabel("$t$")
    ax.set_title("Carta del blow-up $(x,t)$")
    ax.legend(loc="upper left", fontsize=7, ncol=2)

    ax = axs[1]
    for c, col in zip(cs, cols):
        ax.plot(xs, c * xs, color=col, lw=1.1)
    ax.plot(0, 0, "o", color=C_E1, ms=9, label="$\\pi(E)=\\{0\\}$")
    s = np.linspace(-1.2, 1.2, 400)
    ax.plot(s**2, s**3, color=C_CURVA, lw=1.6, label="$C=\\pi(\\tilde C):\\ y^2=x^3$")
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-2.5, 2.5)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title("Plano $(x,y)$: la recta $t=c$ va a la recta $y=cx$")
    ax.legend(loc="lower right", fontsize=7)
    fig.suptitle("$\\pi:(x,t)\\longmapsto(x,\\,xt)$ contrae $E$ al origen", fontsize=10)
    guardar(fig, "divisor_proyeccion.pdf")
