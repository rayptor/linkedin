import numpy as np
import matplotlib.pyplot as plt

def f(x): 
    return 13*x**9 - x * np.exp(x**7) + np.cos(x) + np.sqrt(8)/x**3 + x**2

x = np.linspace(-2, 2, 1000)
y = f(x)

plt.figure(figsize=(15, 9))
plt.plot(x, y, lw=2, color="red", label="f(x)")
plt.scatter(1.23, 0, color="blue", s=100, zorder=5, label=f"Racine x ≈ {1.23}")
plt.axhline(0, color="k", linestyle="--", lw=1)
plt.axvline(1.23, color="blue", linestyle=":", lw=1.5, alpha=0.7)
plt.xticks(np.arange(-2, 2.1, 0.5))  # Original ticks
plt.ylim(-500, 500)  # Original ylim
plt.xlabel("x")
plt.ylabel("y = f(x)")
plt.legend()
plt.grid(color="k", ls=":", lw=0.5, which='major')
plt.minorticks_on()
plt.show()
