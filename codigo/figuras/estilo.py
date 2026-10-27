"""Estilo común de las figuras (matplotlib por defecto, PDF vectorial, etiquetas en español)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SALIDA = Path(__file__).resolve().parent / "figuras_salida"
SALIDA.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42,      # fuentes TrueType embebidas: texto vectorial y seleccionable
    "font.size": 9,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
})

# Colores consistentes con el apéndice
C_CURVA = "tab:red"     # curva / transformada estricta
C_E1 = "tab:blue"
C_E2 = "tab:green"
C_E3 = "black"


def guardar(fig, nombre):
    ruta = SALIDA / nombre
    fig.savefig(ruta, metadata={"CreationDate": None})
    print(f"  figura escrita: {ruta}")
