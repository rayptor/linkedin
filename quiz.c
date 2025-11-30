#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdckdint.h>
#include <time.h>

constexpr int X_MIN = 1;
constexpr int X_MAX = 10;
constexpr int QUESTIONS = 7;
constexpr int TEMPSLIMITE = 7;

constexpr size_t zoneTaille = 512;
constexpr size_t droiteTailleMax = zoneTaille - 80;

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof(typeof(arr[0])))

typedef struct {
    char equation[zoneTaille];
    int solution;
} Equation;

static inline int genererEntierAleatoire(int min, int max) {
    return min + rand() % (max - min + 1);
}

static inline bool genererBooleenAleatoire(void) {
    return (rand() & 1) != 0;
}

static inline void genererPartieDroite(
    int valeurCible,
    char* restrict buffer,
    size_t tailleTampon
) {
    if (buffer == nullptr || tailleTampon == 0) {
        return;
    }

    int valeurDepart = genererEntierAleatoire(X_MIN, X_MAX);
    int valeurActuelle = valeurDepart;
    int decalage = snprintf(buffer, tailleTampon, "%d", valeurDepart);

    if (decalage < 0 || (size_t)decalage >= tailleTampon) {
        return;
    }

    int nombreOperations = genererEntierAleatoire(1, 4);

    for (int i = 0; i < nombreOperations; ++i) {
        char operateur = genererBooleenAleatoire() ? '+' : '-';
        int nombre = genererEntierAleatoire(X_MIN, X_MAX);
        int nouvelleValeur;

        if (operateur == '-') {
            if (ckd_sub(&nouvelleValeur, valeurActuelle, nombre) || nouvelleValeur < 1) {
                operateur = '+';
                if (ckd_add(&nouvelleValeur, valeurActuelle, nombre))
                    break;
            }
        } else {
            if (ckd_add(&nouvelleValeur, valeurActuelle, nombre))
                break;
        }

        int zoneEcrite = snprintf(
            buffer + decalage,
            tailleTampon - (size_t)decalage,
            " %c %d",
            operateur,
            nombre
        );

        if (zoneEcrite < 0 || (size_t)decalage + (size_t)zoneEcrite >= tailleTampon)
            return;

        decalage += zoneEcrite;
        valeurActuelle = nouvelleValeur;
    }

    int difference = valeurCible - valeurActuelle;
    if (difference != 0) {
        char op = (difference > 0) ? '+' : '-';
        int absDiff = (difference > 0) ? difference : -difference;
        
        snprintf(
            buffer + decalage,
            tailleTampon - (size_t)decalage,
            " %c %d",
            op,
            absDiff
        );
    }
}

static void genererPartieGauche(Equation* restrict eq) {
    if (eq == nullptr)
        return;

    int coefficientX = genererEntierAleatoire(X_MIN, X_MAX);
    int solutionX = genererEntierAleatoire(X_MIN, X_MAX);
    int valeurCible;

    if (ckd_mul(&valeurCible, coefficientX, solutionX)) {
        coefficientX = X_MIN;
        solutionX = X_MIN;
        valeurCible = coefficientX * solutionX;
    }

    char expressionDroite[droiteTailleMax];
    genererPartieDroite(valeurCible, expressionDroite, sizeof(expressionDroite));

    int taillePartieDroite = snprintf(nullptr, 0, "%d * X = %s", coefficientX, expressionDroite);
    if (taillePartieDroite < 0 || (size_t)taillePartieDroite >= sizeof(eq->equation)) {
        snprintf(eq->equation, sizeof(eq->equation), "%d * X = 1", coefficientX);
        eq->solution = solutionX;
        return;
    }

    snprintf(
        eq->equation,
        sizeof(eq->equation),
        "%d * X = %s",
        coefficientX,
        expressionDroite
    );
    eq->solution = solutionX;
}

static void genererEquationEquilibree(Equation* restrict eq) {
    if (eq == nullptr)
        return;

    int solutionX = genererEntierAleatoire(X_MIN, X_MAX);
    int nombre = genererEntierAleatoire(X_MIN + 1, X_MAX);
    char expressionGauche[64];
    int valeurCible;

    if (genererBooleenAleatoire()) {
        snprintf(expressionGauche, sizeof(expressionGauche), "X * %d", nombre);
        if (ckd_mul(&valeurCible, solutionX, nombre)) {
            solutionX = X_MIN;
            nombre = X_MIN + 1;
            valeurCible = solutionX * nombre;
        }
    } else {
        if (solutionX % nombre == 0) {
            snprintf(expressionGauche, sizeof(expressionGauche), "X / %d", nombre);
            valeurCible = solutionX / nombre;
        } else {
            snprintf(expressionGauche, sizeof(expressionGauche), "X * %d", nombre);
            if (ckd_mul(&valeurCible, solutionX, nombre)) {
                solutionX = X_MIN;
                nombre = X_MIN + 1;
                valeurCible = solutionX * nombre;
            }
        }
    }

    char expressionDroite[droiteTailleMax];
    genererPartieDroite(valeurCible, expressionDroite, sizeof(expressionDroite));

    int taillePartieGauche = snprintf(nullptr, 0, "%s = %s", expressionGauche, expressionDroite);
    if (taillePartieGauche < 0 || (size_t)taillePartieGauche >= sizeof(eq->equation)) {
        snprintf(eq->equation, sizeof(eq->equation), "X * %d = 1", nombre);
        eq->solution = solutionX;
        return;
    }

    snprintf(
        eq->equation,
        sizeof(eq->equation),
        "%s = %s",
        expressionGauche,
        expressionDroite
    );
    eq->solution = solutionX;
}

#pragma GCC diagnostic pop

static inline double calculerTemps(struct timespec debut, struct timespec fin) {
    return (double)(fin.tv_sec - debut.tv_sec) +
           (double)(fin.tv_nsec - debut.tv_nsec) / 1e9;
}

static void afficherSeparateur(void) {
    for (int i = 0; i < 50; ++i) {
        putchar('-');
    }
    putchar('\n');
}

[[nodiscard]] static bool lireReponse(
    char* restrict buffer,
    size_t taille,
    double* restrict tempsEcoule
) {
    struct timespec debut, fin;
    clock_gettime(CLOCK_MONOTONIC, &debut);

    printf("X = ");
    fflush(stdout);

    if (fgets(buffer, (int)taille, stdin) == nullptr)
        return false;

    clock_gettime(CLOCK_MONOTONIC, &fin);
    *tempsEcoule = calculerTemps(debut, fin);

    size_t len = strlen(buffer);
    if (len > 0 && buffer[len - 1] == '\n')
        buffer[len - 1] = '\0';

    return true;
}

static void test(void) {
    printf("\n");
    afficherSeparateur();
    printf("\t\tGÉNÉRATEUR D'ÉQUATIONS CHRONOMÉTRÉ\n");
    afficherSeparateur();
    printf("-> %d équations linéaires à résoudre\n", QUESTIONS);
    printf("-> %d secondes par équation\n", TEMPSLIMITE);
    printf("-> X est un entier dans [%d;%d]\n", X_MIN, X_MAX);
    afficherSeparateur();

    srand((unsigned int)time(nullptr));
    int score = 0;
    double tempsTotal = 0.0;
    char buffer[zoneTaille];

    for (int i = 1; i <= QUESTIONS; ++i) {
        Equation eq = {0};
        if (genererBooleenAleatoire()) {
            genererEquationEquilibree(&eq);
        } else {
            genererPartieGauche(&eq);
        }

        printf("\nQuestion %d/%d :\n", i, QUESTIONS);
        printf("Équation : %s\n", eq.equation);

        double tempsEcoule;
        if (!lireReponse(buffer, sizeof(buffer), &tempsEcoule)) {
            printf("\nErreur de lecture !\n");
            continue;
        }

        char* pointeurFin;
        long reponse = strtol(buffer, &pointeurFin, 10);
        
        bool estValide = !(pointeurFin == buffer || *pointeurFin != '\0');
        bool estCorrect = estValide && (reponse == eq.solution);

        if (!estValide) {
            printf("Faux ! (%.1f s)\n", tempsEcoule);
            printf("La bonne réponse était : X = %d\n", eq.solution);
        } else if (tempsEcoule > TEMPSLIMITE) {
            if (estCorrect) {
                // CORRECTION: Correspond à la logique Python "Trop lent mais correct"
                printf("Trop lent (%.1f s), mais la réponse est correcte !\n", tempsEcoule);
                ++score;
                tempsTotal += tempsEcoule;
            } else {
                // CORRECTION: Correspond à la logique Python "Trop lent et incorrect"
                printf("Trop lent (%.1f s) !\n", tempsEcoule);
                printf("La bonne réponse était : X = %d\n", eq.solution);
            }
        } else {
            // Temps OK
            if (estCorrect) {
                printf("Correct ! (%.1f s sur %ds)\n", tempsEcoule, TEMPSLIMITE);
                ++score;
                tempsTotal += tempsEcoule;
            } else {
                printf("Faux ! (%.1f s)\n", tempsEcoule);
                printf("La bonne réponse était : X = %d\n", eq.solution);
            }
        }
    }

    printf("\n");
    afficherSeparateur();
    printf("\t\tRÉSULTAT FINAL\n");
    afficherSeparateur();

    if (score == QUESTIONS) {
        printf("Tu as calculé %d bonnes réponses en %.1fs "
               "sur les %d questions.\n",
               score, tempsTotal, QUESTIONS);
    } else {
        printf("Tu n'as calculé que %d bonnes réponses "
               "sur les %d questions.\n",
               score, QUESTIONS);
    }

    const char* message;
    if (score == QUESTIONS) {
        message = "C'est parfait !";
    } else if (score == QUESTIONS - 1) {
        message = "C'est très bien !";
    } else if (score >= QUESTIONS / 2) {
        message = "C'est moyen ! Tu peux mieux faire.";
    } else {
        message = "C'est quoi ça ? Retourne donc t'entraîner !";
    }

    printf("%s\n", message);
}

int main(void) {
    test();
    return EXIT_SUCCESS;
}