from sys import argv
import numpy as np
import matplotlib.pyplot as mpl
from sympy import symbols, re, im, lambdify
from sympy.functions.special.spherical_harmonics import Ynm


def harmonique_spherique(n: int, m: int, c: str) -> None:
    symbole_theta = symbols("theta", real=True)
    symbole_phi = symbols("phi", real=True)
    y_lm = Ynm(n, m, symbole_theta, symbole_phi).expand(func=True)
    y_re, y_im = re(y_lm), im(y_lm)

    resolution = 200
    grille_theta = np.linspace(0, np.pi, resolution)
    grille_phi = np.linspace(0, 2 * np.pi, resolution)
    theta, phi = np.meshgrid(grille_theta, grille_phi)
    magnitude = y_re**2 + y_im**2

    if c == "":
        fonction_lambda = lambdify((symbole_theta, symbole_phi), magnitude, modules=["numpy"])
    elif c == "r":
        fonction_lambda = lambdify((symbole_theta, symbole_phi), y_re, modules=["numpy"])
    elif c == "i":
        fonction_lambda = lambdify((symbole_theta, symbole_phi), y_im, modules=["numpy"])
    else:
        raise ValueError("Seules les lettres 'r' et 'i' sont autorisées !")

    rayon = fonction_lambda(theta, phi)
    rayon = np.where(rayon < 0, 0, rayon)
    rho = rayon if c == "" else np.sqrt(rayon)

    x = rho * np.sin(theta) * np.cos(phi)
    y = rho * np.sin(theta) * np.sin(phi)
    z = rho * np.cos(theta)

    ax = mpl.figure(figsize=(8, 8)).add_subplot(111, projection="3d")
    rayon_norme = np.interp(rayon, (rayon.min(), rayon.max()), (0.1, 0.9))
    couleurs = mpl.cm.rainbow(rayon_norme)
    ax.plot_surface(x, y, z, facecolors=couleurs, rstride=2, cstride=1, lw=1, antialiased=True)

    min_max = np.array([[x.min(), x.max()], [y.min(), y.max()], [z.min(), z.max()]])
    centres = np.mean(min_max, axis=1)
    ext = (np.diff(min_max, axis=1) * 0.5).max()
    ax.set_xlim(centres[0] - ext, centres[0] + ext)
    ax.set_ylim(centres[1] - ext, centres[1] + ext)
    ax.set_zlim(centres[2] - ext, centres[2] + ext)
    ax.set_box_aspect([1,1,1])
    ax.set_axis_off()

    mode = ""
    match (c):
        case "r":
            mode = " (réelle)"
        case "i":
            mode = " (imaginaire)"

    mpl.title(f"Harmonique sphérique $Y_{{{n}}}^{{{m}}}(\\theta, \\phi)${mode}")
    mpl.tight_layout()
    # mpl.savefig(f"hs_{n}_{m}{mode}.jpg")
    mpl.show()

if __name__ == "__main__":
    len_args = len(argv)
    c = ""
    if not (3 <= len_args <= 4):
        print(f"Usage : python {argv[0]} l m 'r ou i'")
        exit(1)

    try:
        n = int(argv[1])
        m = int(argv[2])
        if len_args == 4:
            c = str(argv[3])
    except ValueError:
        print("Error: n and m doivent être des entiers.")
        exit(1)

    if np.fabs(m) > n:
        print(f"Erreur : |m| doit être ≤ n : m = {m}, n = {n}")
        exit(1)

    harmonique_spherique(n, m, c)
