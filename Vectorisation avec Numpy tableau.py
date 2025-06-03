import matplotlib.pyplot as plt

ordres = [100, 200, 400, 800, 1600, 3200]
chronos_cholesky = [0.04, 0.33, 2.52, 19.90, 158.69, 1277.60]
chronos_cholesky_vectorisee = [0.02, 0.06, 0.25, 1.03, 4.17, 20.41]

plt.figure(figsize=(12, 7))
plt.plot(ordres,
         chronos_cholesky,
         marker="o",
         ls="-",
         color="#ff3333",
         markersize=7,
         lw=2,
         label="Cholesky")

plt.plot(ordres,
         chronos_cholesky_vectorisee,
         marker="x",
         ls="--",
         color="#33ff00",
         markersize=7,
         lw=2,
         label="Cholesky vectorisée")

police = {"family": "serif", "color":  "k", "weight": "bold", "size": 12 }
plt.xlabel("Ordre de la matrice (100, 200, 400, 800, 1600, 3200)", fontdict=police)
plt.ylabel("Temps d'éxecution en secondes", fontdict=police)
plt.title("Comparaison entre la version boucle et la version vectorisée", fontdict=police, pad=20)
plt.legend(fontsize=10, frameon=True, shadow=True, fancybox=True)
plt.grid(True, which="both", ls=":", lw=0.5, color="gray", alpha=0.7)
plt.tick_params(axis="both", which="major", labelsize=12, grid_linewidth=1, grid_linestyle="-")
plt.tight_layout()
plt.minorticks_on()
plt.gca().set_aspect("equal", adjustable="box")
plt.show()
