import random
import time

MIN_X = 1
MAX_X = 10
QUESTIONS = 7
TEMPS_LIMITE = 7

def generer_partie_droite(valeur_cible: int) -> str:
    nombre_depart: int = random.randint(MIN_X, MAX_X)
    expression_droite: str = str(nombre_depart)
    valeur_actuelle: int = nombre_depart
    operations_droite: int = random.randint(1, 4)

    for _ in range(operations_droite):
        operateur: str = random.choice(['+', '-'])
        nombre: int = random.randint(MIN_X, MAX_X)
        if operateur == '-' and valeur_actuelle - nombre < 1:
            operateur = '+'

        expression_droite += f" {operateur} {nombre}"
        valeur_actuelle += nombre if operateur == '+' else -nombre

    if (difference := valeur_cible - valeur_actuelle) != 0:
        expression_droite += f" {'+' if difference > 0 else '-'} {abs(difference)}"

    return expression_droite


def generer_partie_gauche() -> tuple[str, int]:
    x_coefficient: int = random.randint(MIN_X, MAX_X)
    x_solution: int = random.randint(MIN_X, MAX_X)
    valeur_cible: int = x_coefficient * x_solution
    expression_droite: str = generer_partie_droite(valeur_cible)
    equation: str = f"{x_coefficient} * X = {expression_droite}"

    return equation, x_solution


def generer_equation_equilibree() -> tuple[str, int]:
    x_solution: int = random.randint(MIN_X, MAX_X)
    nombre: int = random.randint(MIN_X+1, MAX_X)
    
    if random.choice([True, False]):
        expression_gauche: str = f"X * {nombre}"
        valeur_cible: int = x_solution * nombre
    else:
        if x_solution % nombre == 0:
            expression_gauche = f"X / {nombre}"
            valeur_cible = x_solution // nombre
        else:
            expression_gauche = f"X * {nombre}"
            valeur_cible = x_solution * nombre

    expression_droite: str = generer_partie_droite(valeur_cible)
    equation: str = f"{expression_gauche} = {expression_droite}"

    return equation, x_solution


def test() -> None:
    print("\n" + "-" * 50)
    print("\t\tGÉNÉRATEUR D'ÉQUATIONS CHRONOMÉTRÉ")
    print("-" * 50)
    print(f"-> {QUESTIONS} équations linéaires à résoudre")
    print(f"-> {TEMPS_LIMITE} secondes par équation")
    print("-> X est un entier dans [1;10]")
    print("-" * 50)
    
    score: int = 0
    limite: int = TEMPS_LIMITE
    nombre_questions: int = QUESTIONS
    temps_total: float = 0
    
    for i in range(1, nombre_questions + 1):
        generateur = random.choice([generer_equation_equilibree, generer_partie_gauche])
        equation: str
        solution: int
        equation, solution = generateur()
        print(f"\nQuestion {i}/{nombre_questions} :")
        print(f"Équation : {equation}")
        debut = time.time()        
        try:
            reponse_utilisateur: str = input("X = ")
            temps: float = time.time() - debut
            est_correct: bool = False
            if reponse_utilisateur.isdigit():
                reponse: int = int(reponse_utilisateur)
                est_correct = (reponse == solution)
            
            if temps > limite:
                if est_correct:
                    print(f"Trop lent ({temps:.1f} s), mais la réponse est correcte !")
                    score += 1
                    temps_total += temps
                else:
                    print(f"Trop lent ({temps:.1f} s) !", end="")
                    print(f" (La bonne réponse était : X = {solution})")
            elif est_correct:
                print(f"Correct ! ({temps:.1f} s sur {limite}s)")
                score += 1
                temps_total += temps
            else:
                print(f"Faux ! ({temps:.1f} s)", end="")
                print(f" (La bonne réponse était : X = {solution})")
                    
        except KeyboardInterrupt:
            print("\n\nQuiz interrompu !")
            break
        except ValueError:
            print("Erreur : Veuillez entrer un nombre entier valide.")
            continue
        except Exception as e:
            print(f"Erreur inattendue : {e}.")
            continue
    
    print("\n" + "-" * 50)
    print("\t\tRÉSULTAT FINAL")
    print("-" * 50)
    
    if score == nombre_questions:
        print(f"Tu as calculé {score} bonnes réponses en {temps_total:.1f}s "
              f"sur les {nombre_questions} questions.")
    else:
        print(f"Tu n'as calculé que {score} bonnes réponses "
              f"sur les {nombre_questions} questions.")

    match score:
        case n if n == nombre_questions:
            message = "C'est parfait !"
        case n if n == nombre_questions - 1:
            message = "C'est très bien !"
        case n if n >= nombre_questions // 2:
            message = "C'est moyen ! Tu peux mieux faire."
        case _:
            message = "C'est quoi ça ? Retourne donc t'entraîner !"

    print(message)

if __name__ == "__main__":
    test()