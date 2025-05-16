from sys import float_info
from math import fabs, sin, cos, atan, sqrt

def newton_rahpson(
        f: callable,
        fd: callable,
        x0: float,
        itmax: int = 20
    ) -> None:
    
    xold = x0
    xnew = xold
    it = 0
    tol = float_info.epsilon

    while it <= itmax:
        fx = f(xnew)
        dfx = fd(xnew)
        if dfx < tol:
            print("Risque de division par zéro, sortie de la fonction.")
            exit
        xnew = xold - fx / dfx
        delta = fabs(xnew - xold)
        if delta < tol:
            print(f"\n--> Racine trouvée, x = {xnew} à {delta} près.")
            break
        print(f"Itération numéro {it:3d} -> x* = {xnew}")
        xold = xnew
        it += 1


def fonction(x: float) -> float:
    ret = 2*x**5 - x**4 * cos(2*x) + x**2 * atan(x + 1) - sin(x) \
        - 2 * sqrt(x + 1) - 1
    return ret

def derivee(x: float) -> float:
    ret = 10*x**4 + 2*x**4 * sin(2*x) - 4*x**3 * cos(2*x) \
        + x**2 / ((x + 1)**2 + 1) \
        - 1 / sqrt(x + 1) - cos(x) + 2*x * atan(x + 1)
    return ret

if __name__ == "__main__":
    x0 = 2
    newton_rahpson(fonction, derivee, x0)
    #x ≈ 1.06792768912482910207653588...

