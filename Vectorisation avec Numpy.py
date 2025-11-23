import numpy as np
import time as t

def chrono(f: callable) -> float:
    def appel(*args, **kargs) -> None:
        debut = t.perf_counter()
        result = f(*args, **kargs)
        fin = t.perf_counter() - debut
        print(f"\t-> Temps de calcul de la fonction {f.__name__}() : {fin:.2f} s") 
        return result
    return appel


def test_matrice(a: np.ndarray) -> None:
    if a.ndim != 2:
        raise np.linalg.LinAlgError("A n'est pas en 2D !")

    if a.size == 0:
        raise np.linalg.LinAlgError("A est vide !")

    n, m = a.shape
    if n != m:
        raise np.linalg.LinAlgError("A est rectangulaire !")
	
    valp = np.linalg.eigvals(a)
    if np.allclose(a, a.T) == False and np.any(valp) <= 0.0:
        raise np.linalg.LinAlgError("A n'est SDP")


@chrono
def cholesky(a: np.array) -> np.ndarray:
    test_matrice(a)
    L = np.copy(a)
    n = L.shape[0]

    for i in range(n):
        tmp = L[i,i]
        L[i,i] = np.sqrt(tmp)

        # éléments de la diagonale
        for j in range(i+1,n):
            L[j,i] /= L[i,i]

        # éléments hors diagonale
        for k in range(i+1,n):
            L[i,k] = 0
            for j in range(k,n):
                L[j,k] -= L[j,i] * L[k,i]

    return L


@chrono
def cholesky_vectorisee(a: np.ndarray) -> np.ndarray:
    test_matrice(a)
    L = np.copy(a)
    n = L.shape[0]
    
    for i in range(n):
        # éléments de la diagonale
        L[i,i] = np.sqrt(a[i,i] - np.sum(L[i,:i]**2))
        
        # éléments hors diagonale
        for j in range(i+1, n):
            L[j,i] = (a[j,i] - np.sum(L[j,:i] * L[i,:i])) / L[i,i]
    
    return L


def creer_matrice_spd(n: int) -> np.ndarray:
    rng = np.random.default_rng()
    m: np.ndarray = rng.random(size=(n,n), dtype=np.float32)

    a = m @ m.T
    d = np.diag(np.diag(a))
    a += d * 2

    return a


if __name__ == "__main__":
    for n in [100, 200, 400, 800, 1600, 3200]:
        matrice = creer_matrice_spd(n)
        print(f"Pour une matrice d'ordre {n} :")
        cholesky(matrice)
        cholesky_vectorisee(matrice)
