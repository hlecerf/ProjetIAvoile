"""
Génère schema_modelisation.png : le schéma de modélisation du problème du voilier
(plan (x, y), départ, port, zones de vent, paramètres alpha et y, point I).
"""

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

Y_PORT, Y_ZONE = 5000, 4000

# Point de virement illustratif (cas y < 4000)
x_v, y_v = 2500, 2800
x_i = x_v * (Y_PORT - Y_ZONE) / (Y_PORT - y_v)
alpha = math.degrees(math.atan2(x_v, y_v))
alpha2 = math.degrees(math.atan2(x_v, Y_PORT - y_v))

fig, ax = plt.subplots(figsize=(9, 7.5))

# Zones de vent
ax.fill_between([0, 3800], 0, Y_ZONE, color="#E6F1FB")
ax.fill_between([0, 3800], Y_ZONE, Y_PORT + 300, color="#E1F5EE")
ax.plot([0, 3800], [Y_ZONE, Y_ZONE], color="gray", ls="--", lw=1)
ax.text(3700, 2000, "Vent fort\ny < 4000 m", ha="right", va="center", color="#0C447C", fontsize=11)
ax.text(3700, 4550, "Vent modéré\ny ≥ 4000 m", ha="right", va="center", color="#085041", fontsize=11)
ax.text(-120, Y_ZONE, "y = 4000", ha="right", va="center", fontsize=10)

# Trajet
ax.annotate("", xy=(x_v, y_v), xytext=(0, 0),
            arrowprops=dict(arrowstyle="-|>", color="#D85A30", lw=2))
ax.annotate("", xy=(0, Y_PORT), xytext=(x_v, y_v),
            arrowprops=dict(arrowstyle="-|>", color="#D85A30", lw=2))
ax.text(1450, 1300, "Bord 1", color="#993C1D", fontsize=10)
ax.text(1450, 3750, "Bord 2", color="#993C1D", fontsize=10)

# Projections du point de virement
ax.plot([x_v, x_v], [0, y_v], ls=":", color="gray", lw=1)
ax.plot([0, x_v], [y_v, y_v], ls=":", color="gray", lw=1)
ax.text(x_v, -230, "x", ha="center", fontsize=11)
ax.text(-120, y_v, "y", ha="right", va="center", fontsize=11)

# Angle alpha au départ (depuis le nord)
ax.plot([0, 0], [0, 900], color="gray", lw=0.8)
ax.add_patch(Arc((0, 0), 1100, 1100, theta1=90 - alpha, theta2=90, color="black"))
ax.text(170, 720, "α", fontsize=14)

# Angle alpha2 au virement (depuis le nord)
ax.plot([x_v, x_v], [y_v, y_v + 800], color="gray", lw=0.8)
ax.add_patch(Arc((x_v, y_v), 900, 900, theta1=90, theta2=90 + alpha2, color="black"))
ax.text(x_v - 290, y_v + 480, "α₂", fontsize=14)

# Points remarquables
ax.plot(0, 0, "ko")
ax.plot(0, Y_PORT, "ko")
ax.plot(x_v, y_v, "o", color="#D85A30")
ax.plot(x_i, Y_ZONE, "ko", ms=5)
ax.text(-120, -80, "Départ (0, 0)", ha="right", fontsize=10)
ax.text(-120, Y_PORT, "Port (0, 5000)", ha="right", va="center", fontsize=10, weight="bold")
ax.text(x_v + 120, y_v - 60, "Virement (x, y)", fontsize=10, weight="bold")
ax.text(x_i + 80, Y_ZONE + 90, "I", fontsize=11, weight="bold")

# Vent
ax.annotate("", xy=(4500, 4100), xytext=(4500, 4900),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=2))
ax.text(4600, 4500, "Vent", va="center", fontsize=10)

# Encadré des paramètres
texte = (
    "Paramètres (gènes)\n"
    "  α ∈ ]20°, 90°[\n"
    "  y ∈ ]0, 5000[\n\n"
    "Déduits\n"
    "  x = y · tan α\n"
    "  α₂ = arctan(x / (5000 − y))\n\n"
    "Contrainte : α₂ > 20°\n"
    "I : passage de y = 4000\n"
    "Virement : +10 s"
)
ax.text(4050, 3300, texte, va="top", fontsize=10, family="DejaVu Sans",
        bbox=dict(boxstyle="round,pad=0.6", fc="white", ec="gray", lw=0.8))

ax.set_xlim(-1700, 6300)
ax.set_ylim(-500, Y_PORT + 400)
ax.set_aspect("equal")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m), nord")
ax.set_title("Modélisation du problème du voilier")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig("schema_modelisation.png", dpi=150)
print("schema_modelisation.png enregistré")
