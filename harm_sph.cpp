#include <cstdlib>
#include <complex>
#include <cmath>
#include <limits>
#include <print>
#include <boost/math/special_functions/spherical_harmonic.hpp>

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::println(stderr,
            "Usage : {} <n> <m> <theta> <phi>\n"
            "       where: n ≥ 0, |m| ≤ n, theta ∈ [0;π], phi ∈ [0;2π] (en radians)",
            argv[0]
        );
        return EXIT_FAILURE;
    }

    int n = std::atoi(argv[1]);
    int m = std::atoi(argv[2]);
    float theta = std::strtof(argv[3], nullptr);
    float phi   = std::strtof(argv[4], nullptr);

    if (n < 0) {
        std::println(stderr, "Erreur : n doit être ≥ 0");
        return EXIT_FAILURE;
    }

    if (std::abs(m) > n) {
        std::println(stderr, "Erreur : |m| doit être ≤ n");
        return EXIT_FAILURE;
    }

    if (theta < 0.0f || theta > static_cast<float>(M_PI) ||
        phi < 0.0f || phi > static_cast<float>(2 * M_PI)) {
        std::println(stderr, "Erreur : theta doit être dans [0, π] et phi dans [0, 2π]");
        return EXIT_FAILURE;
    }

    std::complex<float> Y = boost::math::spherical_harmonic<float>(n, m, theta, phi);

    std::println("Y_{}^{}(theta={:.6f}, phi={:.6f})", n, m, theta, phi);
    std::println("Partie réelle : {:.6f}", Y.real());
    std::println("Partie imaginaire : {:.6f}", Y.imag());

    return EXIT_SUCCESS;
}
