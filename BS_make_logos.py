# -*- coding: utf-8 -*-
"""Generate clean, professional emblem PNGs for the app:
  sai_net_logo.png  - BrainSafe AI emblem (navy rounded square + molecular-brain motif)
  sssihl_logo.png   - tasteful institutional badge PLACEHOLDER (replace with the official logo)
"""
import os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

NAVY = "#0D2137"; NAVY2 = "#163A5F"; GOLD = "#F0A500"; LIGHT = "#EAF2FA"

# ---------- BrainSafe AI emblem ----------
fig, ax = plt.subplots(figsize=(4, 4), dpi=150); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
ax.add_patch(FancyBboxPatch((0.4, 0.4), 9.2, 9.2, boxstyle="round,pad=0.1,rounding_size=1.4",
                            fc=NAVY, ec=GOLD, lw=3))
# molecular-brain motif: nodes + edges (two hemispheres)
nodes = [(3.1,6.6),(4.2,7.4),(5.0,6.4),(4.1,5.5),(3.0,5.4),(3.6,4.4),
         (6.9,6.6),(5.8,7.4),(6.0,5.5),(7.0,5.4),(6.4,4.4)]
edges = [(0,1),(1,2),(2,3),(3,4),(4,0),(3,5),(4,5),
         (6,7),(7,2),(8,2),(8,9),(9,6),(8,10),(9,10),(5,10)]
for a,b in edges:
    ax.plot([nodes[a][0],nodes[b][0]],[nodes[a][1],nodes[b][1]], color=GOLD, lw=1.6, alpha=0.85, zorder=2)
for i,(x,y) in enumerate(nodes):
    ax.add_patch(Circle((x,y), 0.28, fc=(GOLD if i%3 else LIGHT), ec="white", lw=1.2, zorder=3))
ax.text(5.0, 2.35, "BrainSafe", color="white", ha="center", va="center", fontsize=20,
        fontweight="bold", family="DejaVu Sans")
ax.text(5.0, 1.45, "A I", color=GOLD, ha="center", va="center", fontsize=13,
        fontweight="bold", family="DejaVu Sans")
plt.subplots_adjust(0,0,1,1); plt.savefig("sai_net_logo.png", transparent=True, dpi=150); plt.close()

# ---------- SSSIHL institutional badge (PLACEHOLDER) ----------
fig, ax = plt.subplots(figsize=(4, 4), dpi=150); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
ax.add_patch(Circle((5,5), 4.7, fc=NAVY, ec=GOLD, lw=4))
ax.add_patch(Circle((5,5), 3.9, fc="none", ec=GOLD, lw=1.2, alpha=0.7))
# lamp of knowledge motif
ax.add_patch(plt.Polygon([(5,7.2),(4.4,6.2),(5.6,6.2)], closed=True, fc=GOLD, ec="none"))  # flame
ax.add_patch(FancyBboxPatch((4.2,5.2), 1.6, 0.7, boxstyle="round,pad=0.05,rounding_size=0.3", fc=LIGHT, ec="none"))
ax.text(5, 4.4, "S S S I H L", color=LIGHT, ha="center", va="center", fontsize=12, fontweight="bold")
ax.text(5, 3.4, "Sri Sathya Sai Institute", color=GOLD, ha="center", va="center", fontsize=6.6)
ax.text(5, 2.9, "of Higher Learning", color=GOLD, ha="center", va="center", fontsize=6.6)
plt.subplots_adjust(0,0,1,1); plt.savefig("sssihl_logo.png", transparent=True, dpi=150); plt.close()

print("Generated sai_net_logo.png and sssihl_logo.png (SSSIHL is a placeholder - replace with official logo).")
