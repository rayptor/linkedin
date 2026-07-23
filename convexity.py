import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import implicit_multiplication_application, standard_transformations, parse_expr
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import messagebox

EPSILON = -np.finfo(np.float64).eps

def check_convexity():
    try:
        intervalA, intervalB = float(inputA.get()), float(inputB.get())
        if intervalA >= intervalB:
            messagebox.showerror("Erreur", "A doit être < B !")
            return

        expression = tkEntry.get().replace('^', '**')
        transformations = (standard_transformations + (implicit_multiplication_application,))

        x = sp.Symbol('x', real=True)
        f = parse_expr(expression, local_dict={
            'x': x, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan,
            'exp': sp.exp, 'sqrt': sp.sqrt, 'log': sp.log,
            'pi': sp.pi, 'e': sp.E
        }, transformations=transformations)

        d1f, d2f = sp.simplify(sp.diff(f, x, 1)), sp.simplify(sp.diff(f, x, 2))

        xInflection = []
        try:
            roots = sp.solve(d2f, x)
            for r in roots:
                if r.is_real and intervalA <= float(r) <= intervalB:
                    xInflection.append(float(r))
                elif r.is_number and abs(complex(r).imag) < EPSILON and intervalA <= complex(r).real <= intervalB:
                    xInflection.append(complex(r).real)
        except:
            pass

        xInflection = sorted(set(xInflection))

        testPoints = sorted(set([intervalA] + xInflection + [intervalB]))
        intervals = []
        for i in range(len(testPoints) - 1):
            mid = (testPoints[i] + testPoints[i+1]) / 2
            try:
                v = float(d2f.subs(x, mid))
                intervals.append((testPoints[i], testPoints[i+1], "convexe" if v >= EPSILON else "concave"))
            except:
                intervals.append((testPoints[i], testPoints[i+1], "inconnu"))

        def fmt(expr):
            return str(expr).replace('**', '^')

        textOuput.config(state=tk.NORMAL)
        textOuput.delete('1.0', tk.END)
        textOuput.insert(tk.END, f"f(x) = {fmt(f)}\nf'(x)  = {fmt(d1f)}\nf''(x) = {fmt(d2f)}\n\n")

        convexF = [iv for iv in intervals if iv[2] == "convexe"]
        concaveF = [iv for iv in intervals if iv[2] == "concave"]

        if convexF and not concaveF:
            textOuput.insert(tk.END, "CONVEXE sur tout [A, B]\n\n")
        elif concaveF and not convexF:
            textOuput.insert(tk.END, "CONCAVE sur tout [A, B]\n\n")
        else:
            textOuput.insert(tk.END, "Analyse par sous-intervalles :\n\n")
            for a, b, t in intervals:
                textOuput.insert(tk.END, f"[{'CONVEXE' if t == 'convexe' else 'CONCAVE'}] sur [{a:.4f}, {b:.4f}]\n")
            textOuput.insert(tk.END, "\n")

        if xInflection:
            textOuput.insert(tk.END, "Point d'inflexion :\n")
            for ip in xInflection:
                textOuput.insert(tk.END, f"  ({ip:.4f}, {float(f.subs(x, ip)):.4f})\n")

        textOuput.config(state=tk.DISABLED)

        for widget in frameGfx.winfo_children():
            widget.destroy()

        fig, ax = mpl.subplots(1, 1, figsize=(8, 8))
        ax.spines['left'].set_position(('data', 0))
        ax.spines['bottom'].set_position(('data', 0))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', labelsize=10)
        ax.grid(True, alpha=0.3, linestyle='-', lw=0.5)

        xValue = np.linspace(intervalA, intervalB, 2000)
        function = sp.lambdify(x, f, 'numpy')
        functiond2 = sp.lambdify(x, d2f, 'numpy')

        yValue = np.real(function(xValue))
        if np.isscalar(yValue):
            yValue = np.full_like(xValue, yValue, dtype=float)
        else:
            yValue = yValue.astype(float)
        yValue[np.isinf(yValue)] = np.nan

        valsFd2 = np.real(functiond2(xValue))
        if np.isscalar(valsFd2):
            valsFd2 = np.full_like(xValue, valsFd2, dtype=float)
        else:
            valsFd2 = valsFd2.astype(float)
        valsFd2[np.isinf(valsFd2)] = np.nan

        yMedian, yStandard = np.nanmedian(yValue), np.nanstd(yValue)
        if not np.isnan(yMedian) and not np.isnan(yStandard) and yStandard > 0:
            yValue = np.where(np.abs(yValue - yMedian) < 10 * yStandard, yValue, np.nan)

        i = 0
        while i < len(xValue) - 1:
            if np.isnan(valsFd2[i]):
                i += 1
                continue
            isConvex = valsFd2[i] >= EPSILON
            color = "#00FF0D" if isConvex else "#0077FF"
            j = i + 1
            while j < len(xValue) and not np.isnan(valsFd2[j]) and (valsFd2[j] >= EPSILON) == isConvex:
                j += 1
            xSegments, ySegments = xValue[i:j+1], yValue[i:j+1]
            valid = ~np.isnan(ySegments)
            if np.any(valid):
                ax.plot(xSegments[valid], ySegments[valid], color=color, lw=2)
            i = j

        yInflection = []
        for ip in xInflection:
            try:
                yPlot = float(f.subs(x, ip))
                yInflection.append(yPlot)
                ax.plot(ip, yPlot, 'ro', markersize=8, zorder=5)
                ax.plot([ip, ip], [0, yPlot], 'r--', lw=1, alpha=0.7)
                ax.plot([0, ip], [yPlot, yPlot], 'r--', lw=1, alpha=0.7)
            except:
                pass

        ax.set_title(f'f(x) = {fmt(f)}', fontsize=11)
        ax.set_xlabel('x', fontsize=10)
        ax.set_ylabel('y', fontsize=10)

        xMargin = (intervalB - intervalA) * 0.05
        ax.set_xlim(intervalA - xMargin, intervalB + xMargin)

        yValid = yValue[~np.isnan(yValue)]
        if len(yValid) > 0:
            yMin, yMax = min(np.min(yValid), 0), max(np.max(yValid), 0)
            y_margin = max((yMax - yMin) * 0.05, 0.5)
            ax.set_ylim(yMin - y_margin, yMax + y_margin)

        mpl.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=frameGfx)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    except Exception as e:
        messagebox.showerror("Exception -> ", str(e))


root = tk.Tk()
root.title("Test de Convexité avec SymPy")
root.geometry("1000x900")

frameInput = tk.Frame(root, padx=20, pady=20)
frameInput.pack(fill=tk.X)
frameInput.grid_columnconfigure(0, weight=0)
frameInput.grid_columnconfigure(1, weight=1)
frameInput.grid_columnconfigure(2, weight=0)
frameInput.grid_columnconfigure(3, weight=0)

tk.Label(frameInput, text="A =").grid(row=1, column=0, sticky="e")
inputA = tk.Entry(frameInput, width=5, font=('Courier', 11))
inputA.insert(0, "-2")
inputA.grid(row=1, column=1, sticky="w", padx=(2, 0))

tk.Label(frameInput, text="B =").grid(row=2, column=0, sticky="e", pady=(5,0))
inputB = tk.Entry(frameInput, width=5, font=('Courier', 11))
inputB.insert(0, "2")
inputB.grid(row=2, column=1, sticky="w", padx=(2, 0), pady=(5,0))

tk.Label(frameInput, text="f(x) =", font=('Arial', 12)).grid(row=3, column=0, sticky="e", pady=(15,5))

frameInput.grid_columnconfigure(4, weight=1)
tkEntry = tk.Entry(frameInput, width=70, font=('Courier', 12))
tkEntry.insert(0, "-3*x^3 + 2*x^2 + 2*x - 1")
tkEntry.grid(row=3, column=1, columnspan=4, sticky="ew", pady=(15,5))

tk.Button(frameInput, text="Tester la convexite", command=check_convexity, bg="#00FF0D", fg="black", font=('Arial', 12, 'bold'), padx=20, pady=5).grid(row=4, column=0, columnspan=4, pady=15)

textOuput = tk.Text(frameInput, height=8, width=115, font=('Courier', 13), bg="white", fg="black", wrap=tk.NONE, padx=15, pady=15, relief=tk.SOLID, borderwidth=1)
textOuput.grid(row=5, column=0, columnspan=4, sticky='ew')
textOuput.config(state=tk.DISABLED)

hscroll = tk.Scrollbar(frameInput, orient=tk.HORIZONTAL, command=textOuput.xview)
hscroll.grid(row=6, column=0, columnspan=4, sticky='ew')
textOuput.config(xscrollcommand=hscroll.set)

tk.Label(frameInput, text="Fonctions: sin, cos, tan, exp, sqrt, log, +, -, *, /, ^", font=('Arial', 9), fg="gray").grid(row=7, column=0, columnspan=4, sticky='w', pady=(10, 0))
tk.Frame(root, height=2, bg="gray").pack(fill=tk.X, padx=10)

frameGfx = tk.Frame(root, padx=10, pady=10, width=600, height=600)
frameGfx.pack(fill=tk.BOTH, expand=True)
frameGfx.pack_propagate(False)

tk.Label(frameGfx, text="Préparation du Graphique...", font=('Arial', 14), fg="gray").pack(expand=True)

root.mainloop()