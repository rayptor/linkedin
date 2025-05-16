import sys
import cmath
import math

def fitcubic(a: float, b: float, c: float) -> list[complex]:
    ap = a * a
    qq = (ap - 3.0 * b) / 9.0
    rr = (a * (2.0 * ap - 9.0 * b) + 27.0 * c) / 54.0
    qq3 = qq * qq * qq
    rr2 = rr * rr
    inv3 = 1.0 / 3.0

    if rr2 < qq3:
        theta = cmath.acos(rr / (qq ** 1.5))
        cos1 = cmath.cos(theta * inv3)
        gamma = 2.0 * float(cmath.sqrt(qq).real) * float(cos1.real) + a * inv3
    else:
        sqrt_val = cmath.sqrt(rr2 - qq3)
        aa = -math.copysign((abs(rr) + float(sqrt_val.real)) ** inv3, rr)
        bb = 0.0 if aa == 0.0 else qq / aa
        gamma = -aa - bb + a * inv3

    ee = 0.0
    alpha = a - gamma
    beta = b - alpha * gamma
    e1 = 0.0
    e2 = 0.0
    e3 = c - gamma * beta

    for _ in range(16):
        eee = ee
        eeee = eee
        u1 = alpha - gamma
        u2 = beta - gamma * u1
        q1 = e1
        q2 = e2 - gamma * q1
        q3 = e3 - gamma * q2
        delta3 = 0.0 if u2 == 0.0 else q3 / u2
        delta2 = q2 - u1 * delta3
        delta1 = q1 - delta3

        alpha += delta1
        beta += delta2
        gamma += delta3

        e1 = a - gamma - alpha
        e2 = b - alpha * gamma - beta
        e3 = c - gamma * beta
        ee = e1 * e1 + e2 * e2 + e3 * e3

        if math.isclose(ee, eee) or math.isclose(ee, eeee):
            break

    cc1 = alpha * 0.5
    diskr = cc1 * cc1 - beta

    roots = [0j, 0j, 0j]
    if float(diskr.real) >= 0.0:
        sqrt_d = cmath.sqrt(diskr)
        x0 = complex(-cc1 - float(sqrt_d.real), 0.0)
        x1 = complex(beta / x0.real, 0.0)
        roots[0] = x0
        roots[1] = x1
    else:
        sqrt_d = cmath.sqrt(-diskr)
        roots[0] = complex(-cc1, float(sqrt_d.real))
        roots[1] = complex(-cc1, -float(sqrt_d.real))

    roots[2] = complex(-gamma, 0.0)
    return roots

def format_complex(c: complex) -> str:
    real = round(c.real, 6)
    imag = round(c.imag, 6)
    if imag == 0:
        return f"{real}"
    elif imag < 0:
        return f"{real}{imag}i"
    else:
        return f"{real}+{imag}i"

def main():
    if len(sys.argv) != 4:
        print("Usage : python3 script.py a b c")
        return

    try:
        a = float(sys.argv[1])
        b = float(sys.argv[2])
        c = float(sys.argv[3])
    except ValueError:
        print("Erreur : les arguments doivent être des nombres réels.")
        exit

    try:
        roots = fitcubic(a, b, c)
        if roots[0].real == roots[1].real == roots[2].real:
            print("La racine triple est :")
            print(f"x = {format_complex(roots[0])}")
        else:
            print("Les racines sont :")
            for i, root in enumerate(roots, 1):
                print(f"x{i} = {format_complex(root)}")
    except Exception as e:
        print(f"Une erreur est survenue : {e}")

if __name__ == "__main__":
    main()
