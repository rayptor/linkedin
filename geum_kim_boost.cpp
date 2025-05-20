#include <iostream>
#include <limits>
#include <format>
#include <iomanip>
#include <functional>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <boost/multiprecision/mpfr.hpp>

using namespace boost::multiprecision;
using t_mpfr = mpfr_float_100;

template <typename T> T geum_kim(
        std::function<T(T)> f,
        std::function<T(T)> df,
        T x0,
        std::size_t maxit) {

    T xOld = x0;
    T xNew = x0;
    std::size_t k = 0;
    const T tol = std::numeric_limits<T>::epsilon();
    const T beta = static_cast<T>(2);
    const T sigma = -beta;
    const T sigmaSqr = sigma * sigma;
    const T phi1 = static_cast<T>(11) * beta * beta - static_cast<T>(66) * beta
            + static_cast<T>(136);

    while (k < maxit) {
        T fxn = f(xNew);
        T dfxn = df(xNew);

        if (abs(dfxn) < std::numeric_limits<T>::min()) {
            std::cerr << "Dérivée de f(x) trop petite à l'itération "
                    << k << "risque de division par zéro." << std::endl;
            break;
        }

        T yn = xOld - fxn / dfxn;
        T fyn = f(yn);
        T un = (abs(fxn) > tol) ? (fyn / fxn) : static_cast<T>(0);
        T un2 = un * un;
        T phi2 = static_cast<T>(2) * un * (sigmaSqr - static_cast<T>(2) * sigma - static_cast<T>(9))
                - static_cast<T>(4) * sigma - static_cast<T>(6);

        T kfNum = (static_cast<T>(1) + beta * un + (static_cast<T>(-9) + (static_cast<T>(5) * beta)
        / static_cast<T>(2)) * un2);
        T kfDen = (static_cast<T>(1) + (beta - static_cast<T>(2)) * un + (static_cast<T>(-4) + beta
        / static_cast<T>(2)) * un2);
        if (abs(kfDen) < std::numeric_limits<T>::min())
            std::cerr << "Le dénominator kfDen est trop petit à l'itération "
            << k << "risque de division par zéro." << std::endl;
        T kf = kfNum / kfDen;

        T zn = yn - kf * (fyn / dfxn);
        T fzn = f(zn);
        T vn = (abs(fyn) > tol) ? (fzn / fyn) : static_cast<T>(0);
        T wn = (abs(fxn) > tol) ? (fzn / fxn) : static_cast<T>(0);

        T hfNum = (static_cast<T>(1) + static_cast<T>(2) * un + (static_cast<T>(2) + sigma) * wn);
        T hfDen = (static_cast<T>(1) - vn + sigma * wn);
        if (abs(hfDen) < std::numeric_limits<T>::min())
            std::cerr << "Dénominator hfDen trop petit à l'itération "
            << k << "risque de division par zéro." << std::endl;
        T hf = hfNum / hfDen;

        T sn = zn - hf * (fzn / dfxn);
        T fsn = f(sn);
        T tn = (abs(fzn) > tol) ? (fsn / fzn) : static_cast<T>(0);

        T guw1 = static_cast<T>(6) + static_cast<T>(12) * un + static_cast<T>(2) * un2;
        T guw2 = (static_cast<T>(24) - static_cast<T>(11) * beta) + pow(un, 3U) * phi1
                + static_cast<T>(4) * sigma;
        T guw = (static_cast<T>(-0.5)) * un * wn * (guw1 * guw2) + phi2 * wn * wn;

        T wfNum = (static_cast<T>(1) + static_cast<T>(2) * un + (static_cast<T>(2) + sigma) * vn * wn);
        T wfDen = (static_cast<T>(1) - vn - static_cast<T>(2) * wn - tn + static_cast<T>(2)
                * (static_cast<T>(1) + sigma) * vn * wn) + guw;
        if (abs(wfDen) < std::numeric_limits<T>::min())
            std::cerr << "Dénominator wfDen trop petit à l'itération "
                << k << "risque de division par zéro." << std::endl;
        T wf = wfNum / wfDen;

        xNew = sn - wf * (fsn / dfxn);

        T delta = abs(xNew - xOld);
        if (delta < tol) {
            std::cerr << "Convergence terminée :" << std::endl;
            break;
        }
        
        if constexpr (std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::digits10 > 0) {
            std::cout << std::fixed << std::setprecision(std::numeric_limits<T>::digits10)
                << std::format("Itération numéro{:3d} -> X = {}\n", k,
                    xNew.str(std::numeric_limits<T>::digits10, std::ios_base::fixed));
        } else {
            std::cout << std::fixed << std::setprecision(100)
                << std::format("Itération numéro{:3d} -> X = {}\n", k,
                    xNew.str(100, std::ios_base::fixed));
        }

        xOld = xNew;
        ++k;
    }
    if (k == maxit) {
        std::cerr << "Nombre maximal d'itérations " << maxit << " atteint !" << std::endl;
    }

    return xNew;
}

int main(int argc, char* argv[]) {
    t_mpfr valeurInitiale = t_mpfr("1.5");
    std::size_t iterationsMax = 20;
    int precisionSortie = 0;
    std::string msgDef = "Utilisation de la valeur par défaut : ";
    mpfr_float::default_precision(100);
    
    if (argc > 1) {
        try {
            valeurInitiale = t_mpfr(argv[1]);
        } catch (const std::exception& e) {
            std::cerr << "Valeur pour x0 invalide (" << argv[1] << ") !"
                    << msgDef << " '1.5'. Erreur : " << e.what() << "\n";
            valeurInitiale = t_mpfr("1.5");
        }
    }
    if (argc > 2) {
        try {
            long iterationsMaxUser = std::stol(argv[2]);
            if (iterationsMaxUser <= 0) {
                std::cerr << "Nombre maximum d'itérations négatif (" << argv[2] << ")! Doit être > 0. "
                        << msgDef << iterationsMax << std::endl;
                 iterationsMax = 20;
            } else {
                iterationsMax = static_cast<std::size_t>(iterationsMaxUser);
            }
        } catch (const std::invalid_argument& invArg) {
            std::cerr << "Nombre d'itérations invalide (" << argv[2] << ") ! "
                    << msgDef << iterationsMax << "Erreur : " << invArg.what() << std::endl;
            iterationsMax = 20;
        } catch (const std::out_of_range& outOr) {
            std::cerr << "Nombre d'itérations hors limites (" << argv[2] << ") ! "
                    << msgDef << iterationsMax << "Erreur : " << outOr.what() << std::endl;
            iterationsMax = 20;
        }
    }
    std::cout << "Valeur initiale de x0 = " << valeurInitiale << " avec " << iterationsMax
            << " itérations." << std::endl;

    auto f = [&](t_mpfr x) -> t_mpfr {
        return static_cast<t_mpfr>(13) * pow(x, 9U) - x * exp(pow(x, 7U)) + cos(x) + sqrt(static_cast<t_mpfr>(8))
                / pow(x, 3U) + pow(x, 2U);
    };

    auto df = [&](t_mpfr x) -> t_mpfr {
        return static_cast<t_mpfr>(117) * pow(x, 8U) - exp(pow(x, 7U)) * (static_cast<t_mpfr>(7) * pow(x, 7U)
                + static_cast<t_mpfr>(1)) - sin(x) - (static_cast<t_mpfr>(6) * sqrt(static_cast<t_mpfr>(2)))
                / pow(x, 4U) + static_cast<t_mpfr>(2) * x;
    };

    if constexpr (std::numeric_limits<t_mpfr>::is_specialized && std::numeric_limits<t_mpfr>::max_digits10 > 0) {
        precisionSortie = std::numeric_limits<t_mpfr>::max_digits10;
    }

    t_mpfr resultat = geum_kim<t_mpfr>(f, df, valeurInitiale, iterationsMax);
    std::cout << std::fixed << std::setprecision(precisionSortie)
              << std::format("\n\t -> X = {}\n", resultat.str(precisionSortie, std::ios_base::fixed)) << std::endl;

    return 0;
}
