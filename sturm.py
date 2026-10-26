from sympy import Poly, sturm
from sympy.abc import x


def suite_sturm(polynome:Poly, a: float, b: float) -> int:
    suite = sturm(polynome, x)
    
    def changements_de_signe(valeur):
        signes = []
        for p in suite:
            val = p.subs(x, valeur).evalf()
            if val > 0:
                signes.append(1)
            elif val < 0:
                signes.append(-1)
        
        count = 0
        for i in range(len(signes) - 1):
            if signes[i] != signes[i+1]:
                count += 1
        return count

    return changements_de_signe(a) - changements_de_signe(b)


if __name__ == "__main__":
    a, b = -6, 6
    f = x**3 + 3*x**2 - 6*x - 8
    nombre = suite_sturm(f, a, b)
    print(f"\nNombre total de racines réelles de {f} dans ]{a};{b}[ : {nombre}\n")
