"""
Projet Bases de l'I.A. - Trajet optimal d'un voilier par algorithme génétique.

Modélisation :
- Un individu = (alpha, y)
    alpha : cap du bord 1, en degrés, mesuré depuis le nord vers l'est (bord tribord, symétrie)
    y     : ordonnée du point de virement, en mètres
- Grandeurs déduites :
    x      = y * tan(alpha)
    alpha2 = arctan(x / (Y_PORT - y))   (cap du bord 2, mesuré depuis le nord vers l'ouest)
- Coût = temps total du trajet (s), découpé à y = Y_ZONE, + T_VIREMENT, + pénalité si invalide

Algorithme génétique :
- population initiale valide (choix_y)
- sélection des S % meilleurs
- croisement BLX + mutation gaussienne (sigma décroissant)
- enfants invalides pénalisés (pas de réparation)
- élitisme : 1 individu
- arrêt : pas d'amélioration pendant K_STAGNATION générations, ou G_MAX atteint
"""

import math
import random


# =============================================================================
# 1. CONSTANTES DU PROBLÈME
# =============================================================================

Y_PORT = 5000.0        # ordonnée de l'arrivée (m), arrivée en (0, Y_PORT)
Y_ZONE = 4000.0        # frontière vent fort / vent modéré (m)
ALPHA_MIN = 20.0       # en dessous ou égal : vitesse nulle (degrés)
ALPHA_MAX = 90.0       # au-delà le bateau redescend vers le sud (degrés)
T_VIREMENT = 10.0      # temps perdu au virement de bord (s)


# =============================================================================
# 2. PARAMÈTRES DE L'ALGORITHME GÉNÉTIQUE
# =============================================================================

N = 100                # taille de la population (à tester : 50, 100, 500...)
S = 0.3                # proportion de la population sélectionnée comme parents
P_MUTATION = 0.3       # probabilité qu'un enfant soit muté (à comparer avec 1)

BLX_A = 0.25           # extension du croisement BLX : t tiré dans [-BLX_A, 1 + BLX_A]

SIGMA_ALPHA_0 = 1.0    # écart-type initial de la mutation sur alpha (degrés)
SIGMA_Y_0 = 10.0       # écart-type initial de la mutation sur y (m)
DECROISSANCE_SIGMA = 0.95  # sigma *= DECROISSANCE_SIGMA à chaque génération

NB_ELITES = 1          # nombre de meilleurs individus recopiés tels quels

K_STAGNATION = 20      # arrêt si pas d'amélioration pendant K générations
EPSILON = 1e-6         # amélioration minimale (s) pour compter comme une amélioration
G_MAX = 1000           # garde-fou : nombre maximal de générations

PENALITE = 1e6         # temps de base attribué à une solution invalide (s)
COEF_ECART = 1e3       # pénalité supplémentaire par unité d'écart à la contrainte


# =============================================================================
# 3. FONCTIONS PHYSIQUES
# =============================================================================

def vitesse(angle, y):
    """
    Vitesse du bateau (m/s) selon l'angle par rapport au nord (degrés)
    et la zone de vent (déterminée par y).

    - angle <= ALPHA_MIN        -> 0
    - y < Y_ZONE (vent fort)    -> -(120 - angle)**2 / 1000 + 10
    - sinon (vent modéré)       -> -(120 - angle)**2 / 1250 + 8
    """
    if angle <= ALPHA_MIN:
        return 0.0
    if y < Y_ZONE:
        return -(120.0 - angle) ** 2 / 1000.0 + 10.0
    return -(120.0 - angle) ** 2 / 1250.0 + 8.0


def temps_segment(xa, ya, xb, yb, angle):
    """
    Temps (s) pour parcourir le segment A -> B avec un cap 'angle'.
    Le segment est supposé entièrement dans UNE seule zone de vent :
    la zone est déterminée par le milieu du segment.
    Renvoie math.inf si la vitesse est nulle.
    """
    distance = math.hypot(xb - xa, yb - ya)
    if distance == 0.0:
        return 0.0
    v = vitesse(angle, (ya + yb) / 2.0)
    if v <= 0.0:
        return math.inf
    return distance / v


def choix_y(alpha):
    """
    Intervalle des y valides pour un alpha donné, c'est-à-dire ceux qui garantissent alpha2 > ALPHA_MIN.

    Condition : y * tan(alpha) > (Y_PORT - y) * tan(ALPHA_MIN)
      =>  y_min = Y_PORT * tan(ALPHA_MIN) / (tan(alpha) + tan(ALPHA_MIN))
      =>  y dans ]y_min, Y_PORT[

    Renvoie (y_min, y_max).
    """
    t_a = math.tan(math.radians(alpha))
    t_min = math.tan(math.radians(ALPHA_MIN))
    y_min = Y_PORT * t_min / (t_a + t_min)
    return y_min, Y_PORT


# =============================================================================
# 4. CLASSE INDIVIDU
# =============================================================================

class Individu:
    """Un trajet candidat, défini par (alpha, y)."""

    def __init__(self, alpha, y):
        self.alpha = alpha
        self.y = y
        self.cout = None   # mis en cache après le premier calcul

    # --- Grandeurs déduites ---------------------------------------------------

    def x(self):
        """Abscisse du point de virement : x = y * tan(alpha)."""
        return self.y * math.tan(math.radians(self.alpha))

    def alpha2(self):
        """Cap du bord 2 (degrés, depuis le nord vers l'ouest) : arctan(x / (Y_PORT - y))."""
        return math.degrees(math.atan2(self.x(), Y_PORT - self.y))

    def ecart_contrainte(self):
        """
        Mesure de « combien » l'individu viole les contraintes (0 si valide).
        Sert à pénaliser plus fort les solutions très éloignées du domaine valide.
        """
        ecart = 0.0
        if self.alpha <= ALPHA_MIN:
            ecart += ALPHA_MIN - self.alpha + 1.0
        if self.alpha >= ALPHA_MAX:
            ecart += self.alpha - ALPHA_MAX + 1.0
        if self.y <= 0.0:
            ecart += -self.y / 100.0 + 1.0
        if self.y >= Y_PORT:
            ecart += (self.y - Y_PORT) / 100.0 + 1.0
        if ALPHA_MIN < self.alpha < ALPHA_MAX and 0.0 < self.y < Y_PORT:
            a2 = self.alpha2()
            if a2 <= ALPHA_MIN:
                ecart += ALPHA_MIN - a2 + 1.0
        return ecart

    def est_valide(self):
        """
        True si :
        - ALPHA_MIN < alpha < ALPHA_MAX
        - 0 < y < Y_PORT
        - alpha2 > ALPHA_MIN
        """
        return self.ecart_contrainte() == 0.0

    # --- Coût -----------------------------------------------------------------

    def segments(self):
        """
        Liste des segments du trajet : (xa, ya, xb, yb, angle), chacun dans une seule zone.

          - si y < Y_ZONE :
              bord 1  : (0,0) -> (x,y)                  vent fort,   angle alpha
              bord 2a : (x,y) -> I                      vent fort,   angle alpha2
              bord 2b : I -> (0, Y_PORT)                vent modéré, angle alpha2
              avec I = (x * (Y_PORT - Y_ZONE) / (Y_PORT - y), Y_ZONE)
          - si y >= Y_ZONE :
              bord 1a : (0,0) -> J                      vent fort,   angle alpha
              bord 1b : J -> (x,y)                      vent modéré, angle alpha
              bord 2  : (x,y) -> (0, Y_PORT)            vent modéré, angle alpha2
              avec J = (Y_ZONE * tan(alpha), Y_ZONE)
        """
        x, y = self.x(), self.y
        a1, a2 = self.alpha, self.alpha2()
        if y < Y_ZONE:
            xi = x * (Y_PORT - Y_ZONE) / (Y_PORT - y)
            return [
                (0.0, 0.0, x, y, a1),
                (x, y, xi, Y_ZONE, a2),
                (xi, Y_ZONE, 0.0, Y_PORT, a2),
            ]
        xj = Y_ZONE * math.tan(math.radians(a1))
        return [
            (0.0, 0.0, xj, Y_ZONE, a1),
            (xj, Y_ZONE, x, y, a1),
            (x, y, 0.0, Y_PORT, a2),
        ]

    def FonctionCout(self):
        """
        Temps total du trajet (s).
        Si invalide : PENALITE + COEF_ECART * écart à la contrainte.
        """
        ecart = self.ecart_contrainte()
        if ecart > 0.0:
            return PENALITE + COEF_ECART * ecart
        total = T_VIREMENT
        for seg in self.segments():
            total += temps_segment(*seg)
        return total

    def get_cout(self):
        """Renvoie le coût, en le calculant une seule fois."""
        if self.cout is None:
            self.cout = self.FonctionCout()
        return self.cout

    def __repr__(self):
        return f"Individu(alpha={self.alpha:.3f}°, y={self.y:.1f} m)"


# =============================================================================
# 5. OPÉRATEURS GÉNÉTIQUES
# =============================================================================

def generer_population_initiale(n):
    """
    n individus valides :
    - alpha tiré uniformément dans ]ALPHA_MIN, ALPHA_MAX[
    - y tiré uniformément dans choix_y(alpha)
    """
    population = []
    while len(population) < n:
        alpha = random.uniform(ALPHA_MIN, ALPHA_MAX)
        if alpha <= ALPHA_MIN or alpha >= ALPHA_MAX:
            continue
        y_min, y_max = choix_y(alpha)
        y = random.uniform(y_min, y_max)
        ind = Individu(alpha, y)
        if ind.est_valide():
            population.append(ind)
    return population


def selection(population, s):
    """Renvoie la proportion s des meilleurs individus (coût le plus faible), au moins 2."""
    triee = sorted(population, key=lambda ind: ind.get_cout())
    nb = max(2, int(len(triee) * s))
    return triee[:nb]


def croisement_blx(parent1, parent2):
    """
    Croisement BLX : pour chaque paramètre,
        enfant = t * p1 + (1 - t) * p2,  avec t tiré dans [-BLX_A, 1 + BLX_A]
    (un t indépendant par paramètre). Renvoie (alpha, y).
    """
    t_a = random.uniform(-BLX_A, 1.0 + BLX_A)
    t_y = random.uniform(-BLX_A, 1.0 + BLX_A)
    alpha = t_a * parent1.alpha + (1.0 - t_a) * parent2.alpha
    y = t_y * parent1.y + (1.0 - t_y) * parent2.y
    return alpha, y


def mutation_gaussienne(alpha, y, sigma_alpha, sigma_y):
    """
    Avec probabilité P_MUTATION :
        alpha += N(0, sigma_alpha)
        y     += N(0, sigma_y)
    Renvoie (alpha, y).
    """
    if random.random() < P_MUTATION:
        alpha += random.gauss(0.0, sigma_alpha)
        y += random.gauss(0.0, sigma_y)
    return alpha, y


def FaireUnEnfant(parent1, parent2, sigma_alpha, sigma_y):
    """
    Enfant = croisement BLX des deux parents, puis mutation gaussienne.
    Pas de réparation : un enfant invalide sera pénalisé par FonctionCout.
    """
    alpha, y = croisement_blx(parent1, parent2)
    alpha, y = mutation_gaussienne(alpha, y, sigma_alpha, sigma_y)
    return Individu(alpha, y)


def nouvelle_generation(population, sigma_alpha, sigma_y):
    """
    1. Élitisme : recopier les NB_ELITES meilleurs, sans mutation.
    2. Sélection des parents (selection).
    3. Compléter jusqu'à N individus avec FaireUnEnfant(p1, p2),
       p1 et p2 tirés au hasard (distincts) parmi les sélectionnés.
    """
    parents = selection(population, S)
    nouvelle = parents[:NB_ELITES]
    while len(nouvelle) < N:
        p1, p2 = random.sample(parents, 2)
        nouvelle.append(FaireUnEnfant(p1, p2, sigma_alpha, sigma_y))
    return nouvelle


def ConditionArret(historique_meilleurs, generation):
    """
    True si :
    - generation >= G_MAX, ou
    - le meilleur coût ne s'est pas amélioré (de plus de EPSILON) depuis K_STAGNATION générations.

    historique_meilleurs : liste du meilleur coût de chaque génération.
    """
    if generation >= G_MAX:
        return True
    if len(historique_meilleurs) <= K_STAGNATION:
        return False
    ancien = historique_meilleurs[-1 - K_STAGNATION]
    return ancien - historique_meilleurs[-1] < EPSILON


# =============================================================================
# 6. BOUCLE PRINCIPALE
# =============================================================================

def algorithme_genetique(graine=None):
    """Lance l'AG. Renvoie (meilleur individu, historique des meilleurs coûts, nb de générations)."""
    if graine is not None:
        random.seed(graine)

    population = generer_population_initiale(N)
    sigma_alpha, sigma_y = SIGMA_ALPHA_0, SIGMA_Y_0
    historique_meilleurs = []
    generation = 0

    while True:
        meilleur = min(population, key=lambda ind: ind.get_cout())
        historique_meilleurs.append(meilleur.get_cout())

        if ConditionArret(historique_meilleurs, generation):
            break

        population = nouvelle_generation(population, sigma_alpha, sigma_y)
        sigma_alpha *= DECROISSANCE_SIGMA
        sigma_y *= DECROISSANCE_SIGMA
        generation += 1

    return meilleur, historique_meilleurs, generation


# =============================================================================
# 7. VÉRIFICATIONS
# =============================================================================

def tests_manuels():
    """Vérifications sur des cas calculables à la main, avant de lancer l'AG."""
    assert vitesse(20, 0) == 0.0
    assert abs(vitesse(120, 0) - 10.0) < 1e-12
    assert abs(vitesse(120, 4500) - 8.0) < 1e-12
    assert abs(vitesse(20.0001, 0)) < 1e-3

    # choix_y : au bord de l'intervalle, alpha2 vaut exactement ALPHA_MIN
    y_min, _ = choix_y(45.0)
    assert abs(Individu(45.0, y_min).alpha2() - ALPHA_MIN) < 1e-9

    # Pénalités
    assert not Individu(20.0, 2000).est_valide()            # alpha = 20 : vitesse nulle
    assert not Individu(45.0, y_min - 1).est_valide()        # alpha2 < 20
    assert Individu(45.0, y_min + 1).est_valide()
    assert Individu(45.0, y_min - 1).get_cout() >= PENALITE

    # Cas symétrique à la main : alpha = 45°, y = 2500 -> x = 2500, alpha2 = 45°
    ind = Individu(45.0, 2500.0)
    assert abs(ind.x() - 2500.0) < 1e-9 and abs(ind.alpha2() - 45.0) < 1e-9
    v_fort, v_mod = vitesse(45, 0), vitesse(45, 4500)
    d_total_bord2 = math.hypot(2500, 2500)
    d_zone_mod = d_total_bord2 * (Y_PORT - Y_ZONE) / (Y_PORT - 2500)
    attendu = (math.hypot(2500, 2500) / v_fort
               + (d_total_bord2 - d_zone_mod) / v_fort
               + d_zone_mod / v_mod
               + T_VIREMENT)
    assert abs(ind.get_cout() - attendu) < 1e-9

    # Cas y > 4000 : continuité de la géométrie (les segments se raccordent)
    seg = Individu(30.0, 4500.0).segments()
    assert abs(seg[0][2] - seg[1][0]) < 1e-9 and abs(seg[1][2] - seg[2][0]) < 1e-9
    assert abs(seg[2][2]) < 1e-9 and abs(seg[2][3] - Y_PORT) < 1e-9

    print("Tests manuels : OK")


def verification_grille(pas_alpha=0.05, pas_y=1.0):
    """
    Recherche exhaustive sur une grille, pour vérifier le résultat de l'AG (référence pour la partie III).
    Grille grossière puis affinage autour du meilleur point.
    """
    meilleur = None
    alpha = ALPHA_MIN + 0.5
    while alpha < ALPHA_MAX:
        y = 10.0
        while y < Y_PORT:
            ind = Individu(alpha, y)
            if meilleur is None or ind.get_cout() < meilleur.get_cout():
                meilleur = ind
            y += 10.0
        alpha += 0.5

    a0, y0 = meilleur.alpha, meilleur.y
    alpha = a0 - 1.0
    while alpha <= a0 + 1.0:
        y = y0 - 20.0
        while y <= y0 + 20.0:
            ind = Individu(alpha, y)
            if ind.get_cout() < meilleur.get_cout():
                meilleur = ind
            y += pas_y
        alpha += pas_alpha
    return meilleur


# =============================================================================
# 8. GRAPHIQUES
# =============================================================================

def tracer_resultats(meilleur, historique, fichier_convergence="convergence.png",
                     fichier_trajet="trajet_optimal.png"):
    """Courbe de convergence et trajet optimal (nécessite matplotlib)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(historique, color="#185FA5")
    ax.set_xlabel("Génération")
    ax.set_ylabel("Meilleur temps (s)")
    ax.set_title("Convergence de l'algorithme génétique")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fichier_convergence, dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 7))
    ax.axhspan(0, Y_ZONE, color="#E6F1FB", label="Vent fort")
    ax.axhspan(Y_ZONE, Y_PORT + 200, color="#E1F5EE", label="Vent modéré")
    ax.axhline(Y_ZONE, color="gray", ls="--", lw=0.8)
    xs = [0.0, meilleur.x(), 0.0]
    ys = [0.0, meilleur.y, Y_PORT]
    ax.plot(xs, ys, "-o", color="#D85A30", label="Trajet optimal")
    ax.annotate(f"Virement\n({meilleur.x():.0f}, {meilleur.y:.0f})",
                (meilleur.x(), meilleur.y), xytext=(10, 0), textcoords="offset points")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Trajet optimal : {meilleur.get_cout():.1f} s")
    ax.set_xlim(-200, max(xs) + 1600)
    ax.set_ylim(-200, Y_PORT + 200)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(fichier_trajet, dpi=150)
    plt.close(fig)


# =============================================================================
# 9. LANCEMENT
# =============================================================================

if __name__ == "__main__":
    tests_manuels()

    meilleur, historique, nb_gen = algorithme_genetique(graine=42)
    print(f"Meilleur trajet : {meilleur}")
    print(f"x = {meilleur.x():.1f} m, alpha2 = {meilleur.alpha2():.3f}°")
    print(f"Temps total : {meilleur.get_cout():.3f} s en {nb_gen} générations")

    ref = verification_grille()
    print(f"Référence (recherche sur grille) : {ref}, temps = {ref.get_cout():.3f} s")
    print(f"Écart AG / grille : {meilleur.get_cout() - ref.get_cout():+.3f} s")

    try:
        tracer_resultats(meilleur, historique)
        print("Graphiques enregistrés : convergence.png, trajet_optimal.png")
    except ImportError:
        print("matplotlib non installé : pas de graphiques.")
