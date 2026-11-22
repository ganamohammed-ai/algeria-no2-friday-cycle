#!/usr/bin/env python3
"""Figure 3: Friday effect on the NO2 column (% relative to Sunday-Thursday), three models, with 95% CIs.
Usage: python fig_friday.py --effects friday_cycle.csv --out fig_friday_effect.png"""
import argparse, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
NAMES = {"algiers": "Algiers", "oran": "Oran", "arzew": "Arzew", "constantine": "Constantine", "setif": "Sétif",
         "hassi_messaoud": "Hassi Messaoud", "control_tell": "Control (Tell)", "control_highland": "Control (highland)"}
SER = [("raw", "Unadjusted", "#1f5fa8", "o", -0.22), ("adj", "Weather-adjusted", "#c2591a", "s", 0.0),
       ("ring", "Weather + regional background", "#3a8a5c", "D", 0.22)]
ap = argparse.ArgumentParser(); ap.add_argument("--effects", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
rows = {r["site"]: r for r in csv.DictReader(open(a.effects))}
order = list(NAMES); ypos = {}; y = 0
for i, s in enumerate(order):
    if i == 6: y += 0.7
    ypos[s] = y; y += 1.2
fig, ax = plt.subplots(figsize=(8.4, 5.8), dpi=200)
for tag, lab, col, mk, dy in SER:
    for k, s in enumerate(order):
        r = rows[s]; v, lo, hi = float(r["Fri_%s_pct" % tag]), float(r["Fri_%s_lo" % tag]), float(r["Fri_%s_hi" % tag])
        ax.plot([lo, hi], [ypos[s] + dy] * 2, color=col, lw=1.8, solid_capstyle="round", zorder=2)
        ax.plot([v], [ypos[s] + dy], marker=mk, color=col, ms=6.5, mec="#fcfcfb", mew=1.0, zorder=3, label=lab if k == 0 else None)
ax.axvline(0, color="#2b2b2b", lw=1, zorder=1)
ax.axhline(ypos["hassi_messaoud"] + 0.95, color="#e3e3df", lw=1)
ax.set_yticks([ypos[s] for s in order]); ax.set_yticklabels([NAMES[s] for s in order], fontsize=10, color="#2b2b2b")
ax.invert_yaxis(); ax.grid(axis="x", color="#e3e3df", lw=0.8, zorder=0); ax.set_axisbelow(True)
ax.set_xlabel("Friday vs Sunday–Thursday: change in NO₂ column (%, 95% CI)", fontsize=10, color="#2b2b2b")
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(a.out, facecolor="#fcfcfb")
