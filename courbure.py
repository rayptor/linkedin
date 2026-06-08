# https://docs.sympy.org/latest/modules/diffgeom.html
import sympy as sp
from sympy import sqrt, simplify, diff, Matrix
from sympy.matrices.dense import MutableDenseMatrix

u, v = sp.symbols("u v", real=True)


class FundamentalForms:
    __slots__ = ("E", "F", "G", "L", "M", "N")

    def __init__(
        self,
        E: sp.Expr, F: sp.Expr, G: sp.Expr,
        L: sp.Expr, M: sp.Expr, N: sp.Expr
    ) -> None:
        # First fundamental form: I = E du² + 2F du dv + G dv²
        self.E, self.F, self.G = E, F, G
        # Second fundamental form: II = L du² + 2M du dv + N dv²
        self.L, self.M, self.N = L, M, N


class EvaluatedCurvatures:
    __slots__ = ("point", "K", "H", "k1", "k2")

    def __init__(
        self,
        point: str,
        K: sp.Expr, H: sp.Expr, k1: sp.Expr, k2: sp.Expr
    ) -> None:
        self.point = point
        self.K, self.H, self.k1, self.k2 = K, H, k1, k2


class CurvatureResult:
    __slots__ = ("name", "equation", "forms", "K", "H", "k1", "k2", "notes")

    def __init__(
        self,
        name: str,
        equation: str,
        forms: FundamentalForms,
        K: sp.Expr, H: sp.Expr, k1: sp.Expr, k2: sp.Expr,
        notes: str = ""
    ) -> None:
        
        self.name, self.equation = name, equation
        self.forms = forms
        self.K, self.H, self.k1, self.k2 = K, H, k1, k2
        self.notes = notes


    def evaluate_at(self, substitutions: dict) -> EvaluatedCurvatures:
        def sub(expr: sp.Expr) -> sp.Expr:
            return simplify(expr.subs(substitutions))

        return EvaluatedCurvatures(
            point=", ".join(f"{k} = {val}" for k, val in substitutions.items()),
            K = sub(self.K),
            H = sub(self.H),
            k1 = sub(self.k1),
            k2 = sub(self.k2)
        )


def normal_vector(
        p_u: MutableDenseMatrix,
        p_v: MutableDenseMatrix
) -> MutableDenseMatrix:
    
    cross = p_u.cross(p_v)
    return cross / sqrt(cross.dot(cross))


def compute_fundamental_forms(
    p_u: MutableDenseMatrix,
    p_v: MutableDenseMatrix
) -> FundamentalForms:
    
    E = simplify(p_u.dot(p_u))
    F = simplify(p_u.dot(p_v))
    G = simplify(p_v.dot(p_v))

    n = normal_vector(p_u, p_v)
    p_uu = Matrix([diff(p_u[i], u) for i in range(3)])
    p_uv = Matrix([diff(p_u[i], v) for i in range(3)])
    p_vv = Matrix([diff(p_v[i], v) for i in range(3)])

    L = simplify(n.dot(p_uu))
    M = simplify(n.dot(p_uv))
    N = simplify(n.dot(p_vv))

    return FundamentalForms(E = E, F = F, G = G, L = L, M = M, N = N)


def compute_curvatures(
    forms: FundamentalForms
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    E, F, G = forms.E, forms.F, forms.G
    L, M, N = forms.L, forms.M, forms.N

    det_I = E * G - F**2
    K  = simplify((L * N - M**2) / det_I)
    H  = simplify((E * N - 2 * F * M + G * L) / (2 * det_I))
    delta = sp.sqrt(simplify(H**2 - K))
    k1 = simplify(H + delta)
    k2 = simplify(H - delta)

    return K, H, k1, k2


def surface_analysis(
    name: str,
    equation: str,
    parametrisation: MutableDenseMatrix,
    *,
    notes: str = ""
) -> CurvatureResult:
    
    p_u = Matrix([diff(parametrisation[i], u) for i in range(3)])
    p_v = Matrix([diff(parametrisation[i], v) for i in range(3)])
    forms = compute_fundamental_forms(p_u, p_v)
    K, H, k1, k2 = compute_curvatures(forms)
    return CurvatureResult(
        name = name, equation = equation, forms = forms,
        K = K, H = H, k1 = k1, k2 = k2, notes = notes
    )


_SEPARATOR = "=" * 60

def fmt(expr: sp.Expr) -> str:
    return str(simplify(expr))


def print_curvature_report(
    result: CurvatureResult
) -> None:
    print(f"\n{_SEPARATOR}")
    print(f"{result.name} ({result.equation})")
    print(_SEPARATOR)

    f = result.forms
    print("\nFirst fundamental form (I = E du² + 2F du dv + G dv²) :")
    print(_SEPARATOR)
    print(f"\tE = {fmt(f.E)}")
    print(f"\tF = {fmt(f.F)}")
    print(f"\tG = {fmt(f.G)}")

    print("\nSecond fundamental form (II = L du² + 2M du dv + N dv²) :")
    print(_SEPARATOR)
    print(f"\tL = {fmt(f.L)}")
    print(f"\tM = {fmt(f.M)}")
    print(f"\tN = {fmt(f.N)}")

    print("\nCurvatures")
    print(_SEPARATOR)
    print(f"\tGauss K = {fmt(result.K)}")
    print(f"\tMean H = {fmt(result.H)}")
    print(f"\tPrincipal k1 = {fmt(result.k1)}")
    print(f"\tPrincipal k2 = {fmt(result.k2)}")

    if result.notes:
        print(f"\nNote: {result.notes}")



if __name__ == "__main__":
    R = sp.Symbol("R", positive = True)

    sphere = Matrix([
        R * sp.sin(u) * sp.cos(v),
        R * sp.sin(u) * sp.sin(v),
        R * sp.cos(u)
    ])
    print_curvature_report(surface_analysis(
        name="Sphere",
        equation="x² + y² + z² - R² = 0",
        parametrisation=sphere,
        notes="K = 1/R² and H = -1/R : positive constant curvature.\n")
    )
