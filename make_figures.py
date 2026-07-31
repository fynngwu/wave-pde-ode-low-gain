"""Generate placeholder figures for main.tex (fig2.png - fig5.png).

The placeholders mark where the actual simulation figures go.  They can be
regenerated with the real numerical simulation code later.
"""
import os

import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

FIGURES = [
    ("fig2.png", "Example 1: ODE state comparison"),
    ("fig3.png", "Example 1: displacement surface u(x, t)"),
    ("fig4.png", "Example 2: nilpotent ODE state comparison"),
    ("fig5.png", "Example 2: displacement surface u(x, t)"),
]


def make_placeholder(filename, title):
    fig, ax = plt.subplots(figsize=(3.5, 2.6), dpi=150)
    ax.set_facecolor("#f2f2f2")
    ax.text(
        0.5,
        0.58,
        title,
        ha="center",
        va="center",
        fontsize=10,
        weight="bold",
        transform=ax.transAxes,
    )
    ax.text(
        0.5,
        0.38,
        "placeholder - regenerate via simulation",
        ha="center",
        va="center",
        fontsize=8,
        color="#666666",
        transform=ax.transAxes,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#999999")
    fig.tight_layout()
    path = os.path.join(HERE, filename)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    for filename, title in FIGURES:
        make_placeholder(filename, title)
