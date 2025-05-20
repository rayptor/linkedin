# mypy --ignore-missing-imports geum_kim_mpmath.py
import sys
from mpmath import mp, mpf, sin, cos, sqrt, exp, power, fabs
from typing import Callable, Any

mp.dps = 100

def geum_kim(
        f: Callable[[Any], Any],
        df: Callable[[Any], Any],
        x0: Any,
        maxit: int,
    ) -> Any:

    xOld = mpf(x0)
    xNew = mpf(xOld)
    k = 0
    tol = power(mpf(10), -(mp.dps - 5))
    beta = mpf(2)
    sigma = -beta
    sigmaSqr = power(sigma, 2)
    phi1 = mpf(11) * power(beta,2) - mpf(66) * beta + mpf(136)

    while k < maxit:
        fxn = f(xNew)
        dfxn = df(xNew)
        
        if fabs(dfxn) < mp.eps:
            print(f"Dérivée de f(x) trop petite à l'itération {k}.")
            break
            
        yn = xOld - fxn / dfxn
        fyn = f(yn)
        un = mpf(0) if fabs(fxn) < tol else fyn / fxn            
        un2 = power(un, 2)
        phi2 = mpf(2) * un * (sigmaSqr - mpf(2) * sigma - mpf(9)) - mpf(4) * sigma - mpf(6)

        kfNum = (mpf(1) + beta * un + (mpf(-9) + (mpf(5) * beta) / mpf(2)) * un2)
        kfDen = (mpf(1) + (beta - mpf(2)) * un + (mpf(-4) + beta / mpf(2)) * un2)
        if fabs(kfDen) < mp.eps:
            print(f"Le dénominator kfDen est trop petit à l'itération {k}.")
        kf = kfNum / kfDen

        zn = yn - kf * (fyn / dfxn)
        fzn = f(zn)
        vn = mpf(0) if fabs(fyn) < tol else fzn / fyn
        wn = mpf(0) if fabs(fxn) < tol else fzn / fxn
        
        hfNum = (mpf(1) + mpf(2) * un + (mpf(2) + sigma) * wn)
        hfDen = (mpf(1) - vn + sigma * wn)
        if fabs(hfDen) < mp.eps:
            print(f"Dénominator hfDen trop petit à l'itération {k}.")
            break
        hf = hfNum / hfDen

        sn = zn - hf * (fzn / dfxn)
        fsn = f(sn)
        tn = mpf(0) if fabs(fzn) < tol else fsn / fzn
            
        guw1 = mpf(6) + mpf(12) * un + mpf(2) * un2
        guw2 = (mpf(24) - mpf(11) * beta) + power(un,3) * phi1 + mpf(4) * sigma
        guw = (mpf(-0.5)) * (un * wn * (guw1 * guw2)) + phi2 * power(wn,2)
        
        wfNum = (mpf(1) + mpf(2) * un + (mpf(2) + sigma) * vn * wn)
        wfDen = (mpf(1) - vn - mpf(2) * wn - tn + mpf(2) * (mpf(1) + sigma) * vn * wn) + guw
        if fabs(wfDen) < mp.eps:
            print(f"Dénominator wfDen trop petit à l'itération {k} risque de division par zéro.")
            break
        wf = wfNum / wfDen

        xNew = sn - wf * (fsn / dfxn)

        delta = fabs(xNew - xOld)
        if delta < tol:
            print(f"Convergence terminée :")
            break
        
        print(f"Itération numéro{k:3d} -> X = {xNew}")
        xOld = xNew
        k += 1
    else:
        print(f"Nombre maximal d'itérations '{maxit}' atteint !")

    return xNew


if __name__ == "__main__":
    valeurInitiale = mpf("1.5")
    iterationsMax = int(20)
    iterationsMaxUser = int(0)
    msgDef = "Utilisation de la valeur par défaut : "

    if len(sys.argv) > 1:
        try:
            valeurInitiale = mpf(sys.argv[1])
        except Exception as e:
            print(f"Valeur pour x0 invalide ({sys.argv[1]}) ! {msgDef} '1.5'. Erreur : {e}", file=sys.stderr)
            valeurInitiale = mpf("1.5")
            
    if len(sys.argv) > 2:
        try:
            iterationsMaxUser = int(sys.argv[2])
            if iterationsMaxUser <= 0:
                print(f"Nombre maximum d'itérations négatif ({sys.argv[2]}) ! Doit être > 0. \
                      {msgDef}{iterationsMax}", file=sys.stderr)
                iterationsMax = 20
            else:
                iterationsMax = iterationsMaxUser
        except ValueError as invArg:
            print(f"Nombre d'itérations invalide ({sys.argv[2]}) ! {msgDef}{iterationsMax} \
                  Erreur : {invArg}", file=sys.stderr)
            iterationsMax = 20
        except OverflowError as outOr:
            print(f"Nombre d'itérations hors limites ({sys.argv[2]}) ! {msgDef}{iterationsMax} \
                  Erreur : {outOr}", file=sys.stderr)
            iterationsMax = 20

    print(f"Valeur initiale de x0 = {valeurInitiale} avec {iterationsMax} itérations.")
    
    f = lambda x: mpf(13)*power(x,9) - x * exp(power(x,7)) + cos(x) + sqrt(mpf(8))/power(x,3) + power(x,2)
    df = lambda x: mpf(117)*power(x,8) - exp(power(x,7))*(mpf(7)*power(x,7) + mpf(1)) - sin(x) - \
        (mpf(6)*sqrt(mpf(2)))/power(x,4) + mpf(2)*x

    resultat = geum_kim(f, df, valeurInitiale, iterationsMax)    
    print(f"\n\t -> X = {resultat}\n")
