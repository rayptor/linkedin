// g++ -march=native -std=c++23 -O3 -Wall -pedantic fitcs_strobach.cpp -o fitcs_strobach

#include <array>
#include <cmath>
#include <complex>
#include <format>
#include <expected>
#include <print>
#include <string>
#include <cstdlib>

void fitcs(double a, double b, double c, std::array<std::complex<double>, 3>& x) {
    const auto ap = a * a;
    const auto qq = (ap - 3.0 * b) / 9.0;
    const auto rr = (a * (2.0 * ap - 9.0 * b) + 27.0 * c) / 54.0;
    const auto qq3 = qq * qq * qq;
    const auto rr2 = rr * rr;
    const double inv3 = 1.0 / 3.0;
    double gamma;

    if (rr2 < qq3) {
        const auto theta = std::acos(rr / std::pow(qq, 1.5));
        const auto cos1 = std::cos(theta * inv3);
        gamma = std::fma(2.0 * std::sqrt(qq), cos1, a * inv3);
    } else {
        const auto aa = -std::copysign(std::pow(std::abs(rr) \
                        + std::sqrt(rr2 - qq3), inv3), rr);
        const auto bb = (aa == 0.0) ? 0.0 : qq / aa;
        gamma = -aa - bb + a * inv3;
    }

    double ee = 0.0;
    double alpha = a - gamma;
    double beta = b - alpha * gamma;
    double e1 = 0.0, e2 = 0.0, e3 = c - gamma * beta;

    for (int k = 0; k < 16; ++k) {
        double eee = ee;
        double eeee = eee;
        double u1 = alpha - gamma;
        double u2 = beta - gamma * u1;
        double q1 = e1;
        double q2 = e2 - gamma * q1;
        double q3 = e3 - gamma * q2;
        double delta3 = (u2 == 0.0) ? 0.0 : q3 / u2;
        double delta2 = q2 - u1 * delta3;
        double delta1 = q1 - delta3;

        alpha += delta1;
        beta += delta2;
        gamma += delta3;

        e1 = a - gamma - alpha;
        e2 = b - alpha * gamma - beta;
        e3 = c - gamma * beta;

        ee = std::fma(e1, e1, std::fma(e2, e2, e3 * e3));
        if (ee == 0.0 || ee == eee || ee == eeee)
            break;
    }

    const auto cc1 = alpha * 0.5;
    double diskr = std::fma(cc1, cc1, -beta);

    if (diskr >= 0.0) {
        diskr = std::sqrt(diskr);
        x[0] = {-cc1 - diskr, 0.0};
        x[1] = {beta / x[0].real(), 0.0};
    } else {
        diskr = std::sqrt(-diskr);
        x[0] = {-cc1, diskr};
        x[1] = {-cc1, -diskr};
    }

    x[2] = {-gamma, 0.0};
}

std::string formatx(const std::complex<double>& c) {
    return (c.imag() == 0.0) ? std::format("{}", c.real())
         : std::format("{}{:+}i", c.real(), c.imag());
}

std::expected<std::tuple<double, double, double>, std::string>
arguments(int argc, char* argv[]) {
    if (argc == 4) {
        try {
            double a = std::stod(argv[1]);
            double b = std::stod(argv[2]);
            double c = std::stod(argv[3]);
            return std::tuple{a, b, c};
        } catch (...) {
            return std::unexpected("Erreur de saisie !");
        }
    }
    return std::unexpected("Usage : ./fitcs_strobach a b c");
}

int main(int argc, char* argv[]) {
    auto result = arguments(argc, argv);
    if (!result) {
        std::println("Erreur : {}", result.error());
        std::exit(EXIT_FAILURE);
    }

    auto [a, b, c] = result.value();
    std::array<std::complex<double>, 3> roots;

    try {
        fitcs(a, b, c, roots);

        if (roots[0].real() == roots[1].real() \
        && roots[0].real() == roots[2].real()) {
            std::println("La racine triple est :\nx = {}", formatx(roots[0]));
        } else {
            std::println("Les racines sont :");
            for (std::size_t i = 0; i < roots.size(); ++i)
                std::println("x{} = {}", i + 1, formatx(roots[i]));
        }
    } catch (const std::exception& e) {
        std::println("Une exception est survenue : {}", e.what());
        std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
