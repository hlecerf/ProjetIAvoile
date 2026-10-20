# Projet Bases de l'I.A. : trajet optimal d'un voilier

Recherche, par algorithme génétique, du trajet le plus rapide pour un voilier qui remonte contre le vent de (0, 0) jusqu'au port en (0, 5000), avec un seul virement de bord.

## Contenu

| Fichier | Rôle |
|---|---|
| `Projet IA.pdf` | Consigne du projet |
| `Notes.txt` | Notes de travail et historique des échanges avec le LLM |
| `voilier_ag.py` | Algorithme génétique complet (modèle, opérateurs, tests, vérification) |
| `schema.py` | Génère le schéma de modélisation |
| `schema_modelisation.png` | Schéma du problème dans le plan (x, y) |
| `convergence.png` | Courbe de convergence (meilleur temps par génération) |
| `trajet_optimal.png` | Trajet optimal trouvé |

## Lancer

```bash
pip install matplotlib
python voilier_ag.py   # tests, AG, vérification sur grille, graphiques
python schema.py       # schéma de modélisation
```

## Modélisation

- Individu : (α, y), avec α le cap du premier bord et y l'ordonnée du virement.
- Déduits : x = y·tan α et α₂ = arctan(x / (5000 − y)).
- Coût : temps total, avec les segments découpés à y = 4000, plus 10 s de virement. Les individus invalides (α ≤ 20° ou α₂ ≤ 20°) sont pénalisés.
- Opérateurs : sélection des 30 % meilleurs, croisement BLX, mutation gaussienne à σ décroissant, élitisme sur 1 individu.
- Arrêt : pas d'amélioration pendant 20 générations, ou 1000 générations au maximum.

## Résultat (graine 42)

α = α₂ ≈ 53,64°, virement en (3396, 2500), temps total ≈ 1592,36 s, en une trentaine de générations. Ce résultat est confirmé par une recherche exhaustive sur grille.
