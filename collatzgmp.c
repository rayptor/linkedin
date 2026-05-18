#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

static void collatz_sequence(const mpz_t start, mpz_t summit) {
    mpz_t n, nn;

    mpz_init(nn);
    mpz_init_set(n, start);
    mpz_set(summit, start);

    while (mpz_cmp_ui(n, 1) != 0) {
        if (mpz_even_p(n)) {
            mpz_tdiv_q_2exp(n, n, 1); /* n = n / 2 */
        } else {
            mpz_mul_ui(nn, n, 3); /* n = 3 * n */
            mpz_add_ui(n, nn, 1); /* n = n + 1 */
        }
        if (mpz_cmp(n, summit) > 0) /* n > summit */
            mpz_set(summit, n);
    }

    mpz_clear(n);
    mpz_clear(nn);
}

static void find_altitude_maximale(const mpz_t a, const mpz_t b, mpz_t nMax, mpz_t altMax) {
    mpz_t n, alt;

    mpz_init_set(n, a);
    mpz_init(alt);
    mpz_set(nMax, a);
    mpz_set_ui(altMax, 0);

    while (mpz_cmp(n, b) <= 0) {
        collatz_sequence(n, alt);
        if (mpz_cmp(alt, altMax) > 0) {
            mpz_set(altMax, alt);
            mpz_set(nMax, n);
        }
        mpz_add_ui(n, n, 1); /* n += 1 */
    }
    mpz_clear(alt);
    mpz_clear(n);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <a> <b>  (a < b, strictly positive integers)\n", argv[0]);
        return EXIT_FAILURE;
    }

    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);

    if (mpz_set_str(a, argv[1], 10) != 0 || mpz_set_str(b, argv[2], 10) != 0) {
        fprintf(stderr, "Invalid input(s)!\n");
        mpz_clear(b);
        mpz_clear(a);
        return EXIT_FAILURE;
    }

    if (mpz_sgn(a) <= 0 || mpz_sgn(b) <= 0 || mpz_cmp(a, b) > 0) {
        fprintf(stderr, "A and B must be strictly positive integers such that a < b.\n");
        mpz_clear(b);
        mpz_clear(a);
        return EXIT_FAILURE;
    }

    mpz_t n, altMax;
    mpz_init(n);
    mpz_init(altMax);

    find_altitude_maximale(a, b, n, altMax);

    gmp_printf("\n\tIn [%Zd, %Zd], the maximal altitude '%Zd'\n\tis obtained with the number %Zd.\
        \n\n", a, b, altMax, n);

    mpz_clear(altMax);
    mpz_clear(n);
    mpz_clear(b);
    mpz_clear(a);

    return EXIT_SUCCESS;
}