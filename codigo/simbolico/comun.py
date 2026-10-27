"""Utilidades comunes a los scripts de cálculo simbólico exacto.

Todo el cálculo algebraico de este directorio usa aritmética exacta sobre Q
o sobre cuerpos finitos F_p (sympy); nunca se usan flotantes.
"""
import itertools

from sympy import Poly, groebner, LM, symbols

x, y = symbols("x y")


def encabezado(titulo):
    """Imprime un encabezado legible."""
    print("\n" + "=" * 72)
    print(titulo)
    print("=" * 72)


def comprobar(condicion, mensaje):
    """Comprobación explícita: lanza AssertionError (también con python -O)."""
    if not condicion:
        raise AssertionError("FALLÓ: " + mensaje)
    print("  [OK] " + mensaje)


def exponentes_lider(g, variables, orden="grevlex"):
    """Exponentes del monomio líder de g respecto al orden dado."""
    lm = LM(g, *variables, order=orden)
    return Poly(lm, *variables).monoms()[0]


def dim_cociente(generadores, variables=(x, y), modulus=None, orden="grevlex"):
    """Dimensión sobre el cuerpo base de k[variables]/(generadores).

    Se cuenta el número de monomios estándar de una base de Gröbner.
    Falla (AssertionError) si el cociente no es de dimensión finita.
    """
    G = groebner(list(generadores), *variables, order=orden, modulus=modulus)
    lideres = [exponentes_lider(g, variables, orden) for g in G.exprs]
    cotas = []
    for i in range(len(variables)):
        puras = [e[i] for e in lideres if all(e[j] == 0 for j in range(len(variables)) if j != i)]
        assert puras, "cociente de dimensión infinita"
        cotas.append(min(puras))
    total = 0
    for e in itertools.product(*[range(c) for c in cotas]):
        if not any(all(e[i] >= l[i] for i in range(len(e))) for l in lideres):
            total += 1
    return total


def generadores_potencia_maximal(n, variables=(x, y)):
    """Generadores monomiales de m^n, m = (variables)."""
    gens = []
    for comb in itertools.combinations_with_replacement(range(len(variables)), n):
        m = 1
        for i in comb:
            m *= variables[i]
        gens.append(m)
    return gens
