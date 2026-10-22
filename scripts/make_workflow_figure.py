#!/usr/bin/env python3
"""Study workflow figure for the piRNA breast-cancer diagnostic analysis.

Draws the current design (scripts 06-08):
  7 piRNA cohorts -> label-free preprocessing -> leave-one-cohort-out validation
  (Step 1 feature screening, Step 2 model building & verification)
  -> summary across folds -> functional prediction.

Sample sizes are read from processed_results/ and mRNA_expression/. A cohort whose
file is missing is drawn with a dashed border and "file pending", so re-running
this script after adding the files fills in the numbers.

Usage:   python3 scripts/make_workflow_figure.py
Outputs: results/workflow/Fig_workflow.png / .pdf / .svg
"""
from pathlib import Path
import re
import sys

import logging
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch, FancyBboxPatch,
                                Polygon, Rectangle)
from matplotlib.path import Path as MPath
import pandas as pd

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "workflow"

COHORTS = [  # (file key in processed_results/, display name, specimen)
    ("BRCA1", "TCGA-BRCA", "tissue"),
    ("PRJNA294226", "PRJNA294226", "tissue"),
    ("PRJNA482141", "PRJNA482141", "tissue"),
    ("PRJNA808405", "PRJNA808405", "tissue"),
    ("PRJNA934049", "PRJNA934049", "tissue"),
    ("yyfbatch1", "yyfbatch1", "plasma"),
    ("yyfbatch2", "yyfbatch2", "plasma"),
]
TUMOR = {"Tumor", "Cancer", "cancer", "tumor"}
NORMAL = {"Normal", "Benign", "benign", "normal", "control", "Control", "healthy", "Healthy"}
N_FEATURES = 10          # K_FEATURES in 08_leave_one_cohort_out.R
N_FOLDS = len(COHORTS)


# ------------------------------------------------------------------ data ----
def cohort_file(key):
    for suffix in ("_processed.csv", "_processed_1.csv"):
        f = ROOT / "processed_results" / f"{key}{suffix}"
        if f.exists():
            return f
    return None


def cohort_counts(key):
    f = cohort_file(key)
    if f is None:
        return None
    g = pd.read_csv(f, usecols=["Group"])["Group"].astype(str)
    return int(g.isin(TUMOR).sum()), int(g.isin(NORMAL).sum())


def mrna_info():
    f = ROOT / "mRNA_expression" / "TCGA_BRCA_mRNA.csv"
    p = cohort_file("BRCA1")
    if not f.exists() or p is None:
        return None
    with open(f) as fh:
        genes = len(set(fh.readline().rstrip("\n").split("\t")[1:]))
    mrna_ids = set(pd.read_csv(f, sep="\t", usecols=[0]).iloc[:, 0].astype(str))
    pirna_ids = set(pd.read_csv(p, usecols=[0]).iloc[:, 0].astype(str))
    return len(mrna_ids & pirna_ids), genes


# ----------------------------------------------------------------- style ----
plt.rcParams.update({
    "font.family": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "svg.fonttype": "none",       # keep text editable in the SVG
    "pdf.fonttype": 42,           # embed TrueType so PDF text stays editable
})
INK, INK2, MUTED = "#1f1f1e", "#45443f", "#6b6a64"
ARROW = "#8a897f"
DATA, PREP, S1, S2, FOLD, FUNC = ("#2a78d6", "#6b6a64", "#eb6834", "#1baf7a",
                                  "#eda100", "#e87ba4")


def tint(c, a):
    """Mix colour c with white; a = weight of the colour (0..1)."""
    r, g, b = to_rgb(c)
    return (1 - (1 - r) * a, 1 - (1 - g) * a, 1 - (1 - b) * a)


W, H = 100.0, 126.0   # drawing units; 1 unit = 0.1 inch
fig = plt.figure(figsize=(W / 10, H / 10))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.set_aspect("equal")
ax.axis("off")

TEXTS = []      # (artist, container rect or None) for the layout checks
ARROWS = []     # axis-aligned polylines for the layout checks


def rbox(x, y, w, h, fc, ec, lw=1.0, ls="-", r=1.0, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, ls=ls, zorder=z))
    return (x, y, x + w, y + h)


def T(x, y, s, size=8.5, weight="normal", color=INK, ha="left", va="center",
      style="normal", box=None, z=6):
    t = ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va,
                style=style, zorder=z)
    TEXTS.append((t, box))
    return t


def arrow(points, color=ARROW, lw=1.3, z=4):
    codes = [MPath.MOVETO] + [MPath.LINETO] * (len(points) - 1)
    ax.add_patch(FancyArrowPatch(path=MPath(points, codes), color=color, lw=lw,
                                 arrowstyle="-|>,head_length=5,head_width=2.8",
                                 mutation_scale=1, shrinkA=0, shrinkB=0, zorder=z,
                                 joinstyle="miter", capstyle="butt"))
    ARROWS.append(points)


def panel(top, h, color, title, fill=0.06, ec=None):
    box = rbox(2, top - h, 96, h, tint(color, fill), ec or color, lw=1.2, r=1.4)
    T(4, top - 2.4, title, 12, "bold", box=box)
    return box


def chip(x, y, w, h, color, lead, details, lead_y=1.9, pitch=1.65, ls="-", center=False):
    box = rbox(x, y, w, h, "white", color, lw=0.9, ls=ls, r=0.8, z=2)
    y_lead = y + h / 2 + pitch * len(details) / 2 + 0.1 if center else y + h - lead_y
    T(x + 1.2, y_lead, lead, 9.5, "bold", box=box)
    for j, d in enumerate(details, 1):
        T(x + 1.2, y_lead - pitch * j, d, 8.5, color=INK2, box=box)
    return box


# ---------------------------------------------------------------- glyphs ----
def cylinder(cx, cy, color, w=2.2, h=2.4):
    body = tint(color, 0.18)
    ax.add_patch(Ellipse((cx, cy - h / 2), w, 0.75, fc=body, ec=color, lw=0.9, zorder=5))
    ax.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h, fc=body, ec="none", zorder=5.1))
    for sx in (cx - w / 2, cx + w / 2):
        ax.add_line(Line2D([sx, sx], [cy - h / 2, cy + h / 2], color=color, lw=0.9, zorder=5.2))
    ax.add_patch(Ellipse((cx, cy + h / 2), w, 0.75, fc=tint(color, 0.4), ec=color, lw=0.9,
                         zorder=5.3))


def tube(cx, cy, color, w=1.15, h=2.9):
    y0 = cy - h / 2
    ax.add_patch(FancyBboxPatch((cx - w / 2, y0), w, h * 0.78,
                                boxstyle=f"round,pad=0,rounding_size={w / 2 - 0.01}",
                                fc="white", ec=color, lw=0.9, zorder=5))
    ax.add_patch(FancyBboxPatch((cx - w / 2 + 0.12, y0 + 0.12), w - 0.24, h * 0.36,
                                boxstyle=f"round,pad=0,rounding_size={w / 2 - 0.13}",
                                fc=tint(color, 0.45), ec="none", zorder=5.1))
    ax.add_patch(Rectangle((cx - w * 0.62, y0 + h * 0.74), w * 1.24, h * 0.26,
                           fc=color, ec=color, lw=0.6, zorder=5.2))


def forest(x, y, w, h, color):
    to_x = lambda auc: x + w * (auc - 0.3) / 0.7
    ax.add_line(Line2D([to_x(0.5)] * 2, [y, y + h], color=MUTED, lw=0.8,
                       ls=(0, (2, 1.6)), zorder=5))
    for i, (est, lo, hi) in enumerate([(0.74, 0.57, 0.89), (0.62, 0.47, 0.77),
                                       (0.80, 0.66, 0.92)]):
        yy = y + h * (0.82 - 0.32 * i)
        ax.add_line(Line2D([to_x(lo), to_x(hi)], [yy, yy], color=color, lw=1.6, zorder=5))
        d = 0.65
        ax.add_patch(Polygon([(to_x(est) - d, yy), (to_x(est), yy + d),
                              (to_x(est) + d, yy), (to_x(est), yy - d)],
                             fc=color, ec="white", lw=0.6, zorder=5.1))


def network(cx, cy, s, color):
    pi = (cx - s, cy)
    genes = [(cx, cy + 0.62 * s), (cx, cy), (cx, cy - 0.62 * s)]
    pws = [(cx + s, cy + 0.36 * s), (cx + s, cy - 0.36 * s)]
    for g in genes:
        ax.add_line(Line2D([pi[0], g[0]], [pi[1], g[1]], color=ARROW, lw=0.8, zorder=5))
    for g, p in [(0, 0), (1, 0), (1, 1), (2, 1)]:
        ax.add_line(Line2D([genes[g][0], pws[p][0]], [genes[g][1], pws[p][1]],
                           color=tint(color, 0.8), lw=0.8, zorder=5))
    d = 0.75
    ax.add_patch(Polygon([(pi[0] - d, pi[1]), (pi[0], pi[1] + d), (pi[0] + d, pi[1]),
                          (pi[0], pi[1] - d)], fc=INK, ec="white", lw=0.6, zorder=5.2))
    for g in genes:
        ax.add_patch(Circle(g, 0.5, fc=color, ec="white", lw=0.6, zorder=5.2))
    for p in pws:
        ax.add_patch(Rectangle((p[0] - 0.6, p[1] - 0.6), 1.2, 1.2, fc=tint(color, 0.45),
                               ec=color, lw=0.8, zorder=5.2))


# ============================================================ A. data ======
counts = {key: cohort_counts(key) for key, _, _ in COHORTS}
mrna = mrna_info()

top = H - 2.0
hA = 19.5
A = panel(top, hA, DATA, "Data collection")
tile_w, tile_h, gap, ggap = 11.7, 9.4, 1.0, 3.0
row_w = 7 * tile_w + 5 * gap + ggap
xs, x = [], 2 + (96 - row_w) / 2
for i in range(7):
    xs.append(x)
    x += tile_w + (ggap if i == 4 else gap)
tile_y = top - hA + 1.8
label_y = tile_y + tile_h + 1.95
cylinder(xs[0] + 1.1, label_y, DATA)
T(xs[0] + 2.9, label_y, "5 public tissue cohorts", 9.5, "bold", box=A)
tube(xs[5] + 0.7, label_y, DATA)
T(xs[5] + 2.0, label_y, "2 in-house plasma cohorts", 9.5, "bold", box=A)
for (key, name, spec), tx in zip(COHORTS, xs):
    cnt = counts[key]
    pending = cnt is None
    box = rbox(tx, tile_y, tile_w, tile_h, "white", MUTED if pending else DATA, lw=0.9,
               ls=(0, (2.6, 1.8)) if pending else "-", r=0.8, z=2)
    cx, t_top = tx + tile_w / 2, tile_y + tile_h
    T(cx, t_top - 1.9, name, 9, "bold", ha="center", box=box,
      color=INK2 if pending else INK)
    T(cx, t_top - 3.8, spec, 8.5, color=INK2, ha="center", box=box)
    if pending:
        T(cx, t_top - 6.2, "file pending", 8.5, color=MUTED, style="italic", ha="center",
          box=box)
    else:
        T(cx, t_top - 5.75, f"n = {sum(cnt)}", 8.5, ha="center", box=box)
        T(cx, t_top - 7.5, f"{cnt[0]} T / {cnt[1]} N", 8.5, ha="center", box=box)

# ===================================================== B. preprocessing ====
arrow([(50, A[1]), (50, A[1] - 3.5)])
top = A[1] - 3.5
hB = 12.0
B = panel(top, hB, PREP, "Preprocessing", fill=0.05)
steps = [("Harmonise labels", ["Tumor vs Normal"]),
         ("Shared piRNAs", ["present in every cohort"]),
         ("log2(TPM + 1)", ["transformation"]),
         ("Per-cohort z-score", ["label-free, within each cohort"])]
cw, cg = (92 - 3 * 3.4) / 4, 3.4
cy0, ch = top - hB + 1.6, 6.2
for k, (lead, det) in enumerate(steps):
    cx0 = 4 + k * (cw + cg)
    chip(cx0, cy0, cw, ch, PREP, lead, det, pitch=1.85, center=True)
    if k < 3:
        arrow([(cx0 + cw + 0.3, cy0 + ch / 2), (cx0 + cw + cg - 0.3, cy0 + ch / 2)])

# ================================================================ C. LOCO ===
arrow([(50, B[1]), (50, B[1] - 3.5)])
top = B[1] - 3.5
hC = 47.0
C = rbox(2, top - hC, 96, hC, "#f7f7f5", "#a9a89f", lw=1.2, r=1.4)
T(4, top - 2.4, "Leave-one-cohort-out validation", 12, "bold", box=C)
T(4, top - 4.6, f"{N_FOLDS} folds · in each fold one cohort is held out as the independent "
  f"validation set; the other {N_FOLDS - 1} form the discovery set", 9, color=INK2, box=C)
ctop, cbot = top - 9.6, top - hC + 6.5
SW, S1X, S2X = 28.3, 35.3, 66.9

# fold grid: rows = cohorts, columns = folds; fold j holds out cohort j
cell, lab_right = 2.7, 12.0
gx0 = lab_right + 0.6
block_h = 2.6 + 7 * cell + 4.6
g_top = ctop - (ctop - cbot - block_h) / 2 - 2.6
T(lab_right, g_top + 1.25, "Fold", 8.5, "bold", ha="right", box=C)
for j in range(7):
    T(gx0 + (j + 0.5) * cell, g_top + 1.25, str(j + 1), 8.5, ha="center", box=C)
for i, (_, name, _) in enumerate(COHORTS):
    T(lab_right, g_top - (i + 0.5) * cell, name, 8.5, ha="right", box=C)
    for j in range(7):
        ax.add_patch(Rectangle((gx0 + j * cell + 0.14, g_top - (i + 1) * cell + 0.14),
                               cell - 0.28, cell - 0.28, ec="none", zorder=3,
                               fc=S2 if i == j else tint(S1, 0.32)))
g_bot = g_top - 7 * cell
for ly, fc, lab in [(g_bot - 1.7, tint(S1, 0.32), f"Discovery ({N_FOLDS - 1} cohorts)"),
                    (g_bot - 3.7, S2, "Held-out (1 cohort)")]:
    ax.add_patch(Rectangle((gx0, ly - 0.75), 1.5, 1.5, fc=fc, ec="none", zorder=3))
    T(gx0 + 2.1, ly, lab, 8.5, box=C)


def step_box(x, color, label, tagline, title, items):
    box = rbox(x, cbot, SW, ctop - cbot, tint(color, 0.09), color, lw=1.3, r=1.0, z=2)
    pill = rbox(x + 1.3, ctop - 3.1, 6.0, 1.9, "white", color, lw=1.0, r=0.95, z=3)
    T(x + 4.3, ctop - 2.15, label, 8.5, "bold", ha="center", box=pill)
    T(x + 8.1, ctop - 2.15, tagline, 8.5, color=INK2, style="italic", box=box)
    T(x + 1.3, ctop - 4.9, title, 10.5, "bold", box=box)
    yy, lead_ys = ctop - 7.6, []
    for k, (lead, details) in enumerate(items, 1):
        ax.add_patch(Circle((x + 2.3, yy), 0.95, fc="white", ec=color, lw=1.0, zorder=4))
        T(x + 2.3, yy - 0.05, str(k), 8, "bold", ha="center")
        T(x + 4.0, yy, lead, 9, "bold", box=box)
        for j, d in enumerate(details, 1):
            T(x + 4.0, yy - 1.55 * j, d, 8.5, color=INK2, box=box)
        lead_ys.append(yy)
        yy -= 1.55 * (len(details) + 1) + 1.15
    return box, lead_ys


step_box(S1X, S1, "Step 1", "discovery cohorts only", "Feature screening", [
    ("Expression filter", ["detected in ≥ 50% of samples", "of every discovery cohort"]),
    ("Differential expression", ["Welch t-test in each cohort,",
                                 "Tumor vs Normal → signed z"]),
    ("Direction consistency", ["same Tumor/Normal direction",
                               f"in all {N_FOLDS - 1} discovery cohorts"]),
    ("Meta-ranking", ["Stouffer combined |z|", f"→ top {N_FEATURES} candidate piRNAs"]),
])
_, s2_leads = step_box(S2X, S2, "Step 2", "locked before validation",
                       "Model building & verification", [
    ("Logistic regression", ["L2-penalised, class-balanced,",
                             "pooled discovery cohorts"]),
    ("Internal estimate", ["5-fold CV in discovery data,",
                           "Step 1 re-run inside each fold"]),
    ("Lock model + threshold", ["Youden's J on discovery data"]),
    ("Independent validation", ["held-out cohort: AUC (DeLong",
                                "95% CI), sensitivity, specificity"]),
])

# discovery cohorts -> Step 1
g_right = gx0 + 7 * cell
arrow([(g_right + 0.4, g_top - 3.5 * cell), (S1X, g_top - 3.5 * cell)], color=S1, lw=1.5)
# Step 1 -> Step 2
arrow([(S1X + SW, (ctop + cbot) / 2), (S2X, (ctop + cbot) / 2)])
# held-out cohort -> Step 2 item 4 only, routed under Step 1
y_lane, x_v, x_v2 = cbot - 3.0, S1X - 1.9, S2X - 1.65
arrow([(g_right + 0.4, g_top - 6.5 * cell), (x_v, g_top - 6.5 * cell), (x_v, y_lane),
       (x_v2, y_lane), (x_v2, s2_leads[3]), (S2X, s2_leads[3])], color=S2, lw=1.5)
T((x_v + x_v2) / 2, y_lane - 1.35, "held-out cohort → independent validation only", 8.5,
  color=INK2, style="italic", ha="center", box=C)
# loop back to the grid: next fold
y_loop, gx_mid = ctop + 1.7, gx0 + 3.5 * cell
arrow([(S2X + SW / 2, ctop), (S2X + SW / 2, y_loop), (gx_mid, y_loop),
       (gx_mid, g_top + 2.45)], color=MUTED, lw=1.2)
T((gx_mid + S2X + SW / 2) / 2, y_loop + 1.15,
  "next fold — repeat until every cohort has been held out once", 8.5, color=INK2,
  style="italic", ha="center", box=C)

# ===================================================== D. across folds =====
arrow([(50, C[1]), (50, C[1] - 3.5)])
top = C[1] - 3.5
hD = 14.0
D = panel(top, hD, FOLD, f"Across the {N_FOLDS} folds", fill=0.07)
dw, dg = (92 - 2 * 2.6) / 3, 2.6
dy, dh = top - hD + 1.7, 8.2
chip(4, dy, dw, dh, FOLD, "Performance",
     ["discovery-CV vs held-out AUC", "per fold (ROC, forest plot)"], center=True)
forest(4 + dw - 7.6, dy + 1.3, 6.6, dh - 2.6, FOLD)
chip(4 + dw + dg, dy, dw, dh, FOLD, "Selection stability",
     ["how often each piRNA is", "selected across the folds"], center=True)
chip(4 + 2 * (dw + dg), dy, dw, dh, FOLD, "Final signature",
     ["Step 1 and model fit re-run", f"on all {N_FOLDS} cohorts → final model"],
     center=True)

# ================================================== E. functional ==========
arrow([(50, D[1]), (50, D[1] - 3.5)])
top = D[1] - 3.5
hE = 15.5
E = panel(top, hE, FUNC, "Functional prediction of the signature piRNAs", fill=0.07)
ew, eg = (92 - 3 * 3.4) / 4, 3.4
ey, eh = top - hE + 1.7, 9.4
m_txt = (["n = %d samples matched" % mrna[0], "to TCGA piRNA data", f"{mrna[1]:,} genes"]
         if mrna else ["mRNA file pending"])
chip(4, ey, ew, eh, DATA, "TCGA-BRCA mRNA-seq", m_txt)
chip(4 + (ew + eg), ey, ew, eh, FUNC, "piRNA–mRNA correlation",
     ["Pearson r for every gene", "detected in ≥ 70% of samples",
      "FDR < 0.05, |r| > 0.3"])
chip(4 + 2 * (ew + eg), ey, ew, eh, FUNC, "Pathway enrichment",
     ["GO (BP, CC, MF), KEGG,", "Reactome over-representation",
      "GSEA on genes ranked by r"])
chip(4 + 3 * (ew + eg), ey, ew, eh, FUNC, "Functional network", ["piRNA – gene – pathway"])
network(4 + 3 * (ew + eg) + ew / 2, ey + 2.75, 2.8, FUNC)
for k in range(3):
    x0 = 4 + k * (ew + eg) + ew
    arrow([(x0 + 0.3, ey + eh / 2), (x0 + eg - 0.3, ey + eh / 2)])


# ============================================================== checks =====
def layout_problems():
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inv = ax.transData.inverted()

    def bb(t):
        e = t.get_window_extent(rend)
        (x0, y0), (x1, y1) = inv.transform([[e.x0, e.y0], [e.x1, e.y1]])
        return x0, y0, x1, y1

    out, boxes = [], [(bb(t), t.get_text(), c) for t, c in TEXTS]
    for (x0, y0, x1, y1), s, c in boxes:          # text inside its container
        if c and (x0 < c[0] + 0.3 or x1 > c[2] - 0.3 or y0 < c[1] + 0.15 or y1 > c[3] - 0.15):
            out.append(f"overflow: {s!r}")
    for i in range(len(boxes)):                    # text-text overlaps
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][0], boxes[j][0]
            if min(a[2], b[2]) - max(a[0], b[0]) > 0.05 and min(a[3], b[3]) - max(a[1], b[1]) > 0.05:
                out.append(f"overlap: {boxes[i][1]!r} / {boxes[j][1]!r}")
    for pts in ARROWS:                             # arrows crossing text
        for (xa, ya), (xb, yb) in zip(pts, pts[1:]):
            for (x0, y0, x1, y1), s, _ in boxes:
                if (min(xa, xb) < x1 - 0.05 and max(xa, xb) > x0 + 0.05 and
                        min(ya, yb) < y1 - 0.05 and max(ya, yb) > y0 + 0.05):
                    out.append(f"arrow crosses text: {s!r}")
    return out


problems = layout_problems()
for p in problems:
    print("LAYOUT:", p)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUT / f"Fig_workflow.{ext}", dpi=300, facecolor="white")
# make the SVG font portable (Liberation Sans is metric-compatible with Arial)
svg = OUT / "Fig_workflow.svg"
svg.write_text(re.sub(r"font-family:[^;\"]*",
                      "font-family: Arial, Helvetica, 'Liberation Sans', sans-serif",
                      svg.read_text()))
print("Saved:", *(OUT / f"Fig_workflow.{e}" for e in ("png", "pdf", "svg")), sep="\n  ")
print("Cohorts with data:", ", ".join(n for (k, n, _) in COHORTS if counts[k]) or "none")
sys.exit(1 if problems else 0)
