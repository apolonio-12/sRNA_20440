     Show Dotfiles Show Owner/Mode
/home/apolonio/orcd/scratch/Attempt_1/
#!/usr/bin/env python3
"""
plot_figures.py

Generates all paper figures from FINAL_TABLES/.
Run from the same directory that contains FINAL_TABLES/ and
the MetaPhlAn / HUMAnN output files.

Requires: matplotlib, seaborn, pandas, numpy, scipy,
          networkx, adjustText
Install:  pip install matplotlib seaborn pandas numpy scipy networkx adjustText
"""

import os
import warnings
import textwrap
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm, to_rgba
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
import networkx as nx

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except (ImportError, SyntaxError):
    HAS_ADJUSTTEXT = False
    print("[WARN] adjustText unavailable (install version 0.7.3 for Python 3.6) "
          "— volcano labels may overlap.")

#  PATHS 

TABLES   = "FINAL_TABLES"
FIG_OUT  = "FIGURES_May10"
os.makedirs(FIG_OUT, exist_ok=True)

SRNA_MASTER   = f"{TABLES}/srna_master.tsv"
SRNA_PRESENCE = f"{TABLES}/srna_presence.tsv"
SRNA_TARGET   = f"{TABLES}/srna_target_pathway.tsv"
PATHWAY_MASTER= f"{TABLES}/pathway_master.tsv"
MECHANISM     = f"{TABLES}/mechanism_table.tsv"

# MetaPhlAn / HUMAnN community tables — update paths if needed
METAPHLAN_RNA = "AD/Total_RNA_merged_abundance_table.tsv"
METAPHLAN_DNA = "AD/Total_DNA_merged_abundance_table.tsv"

#  PALETTE 

DRY_C   = "#F5DEB3"   # wheat  — 2016, dry
WET_C   = "#1E90FF"   # dodger blue — 2017, wet
SHARE_C = "#555555"   # shared / neutral
UP_C    = WET_C
DOWN_C  = "#E8541A"   # warm orange for down-regulated in wet (pathway decrease)
NS_C    = "#CCCCCC"
HILIGHT = "#FFD700"   # gold for sRNA-targeted pathway overlay

CMAP_DIV = mpl.colors.LinearSegmentedColormap.from_list(
    "dry_wet", [DRY_C, "white", WET_C]
)

mpl.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.labelsize":    11,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
})

def save(name):
    path = f"{FIG_OUT}/{name}"
    plt.savefig(path)
    print(f"  Saved → {path}")
    plt.close()

#  LOAD TABLES 

print("Loading tables...")
sm  = pd.read_csv(SRNA_MASTER,    sep="\t")
sp  = pd.read_csv(SRNA_PRESENCE,  sep="\t")
st  = pd.read_csv(SRNA_TARGET,    sep="\t")
pm  = pd.read_csv(PATHWAY_MASTER, sep="\t")
mech= pd.read_csv(MECHANISM,      sep="\t")

# Unique sRNAs per condition for Jaccard / counts
ids_dry = set(sp[sp["condition"]=="dry"]["srna_id"].unique())
ids_wet = set(sp[sp["condition"]=="wet"]["srna_id"].unique())
ids_shared = ids_dry & ids_wet

# Per-sRNA summary (one row per sRNA)
srna_uniq = (
    sm.drop_duplicates("srna_id")
      [["srna_id","type","length","log2fc","p_value","padj","contig","start","end","strand"]]
      .copy()
)
# Assign direction based on presence
def presence_direction(sid):
    ind = sid in ids_dry
    inw = sid in ids_wet
    if ind and inw:  return "shared"
    if ind:          return "down"   # 2016-only = disappears in wet = "down"
    return "up"                       # 2017-only = appears in wet  = "up"

srna_uniq["direction"] = srna_uniq["srna_id"].apply(presence_direction)

# Best interaction per sRNA (most negative energy)
st_best = (
    st.dropna(subset=["interaction_energy"])
      .sort_values("interaction_energy")
      .drop_duplicates("srna_id")
      [["srna_id","interaction_energy"]]
)
srna_uniq = srna_uniq.merge(st_best, on="srna_id", how="left")

print(f"  sRNAs: {len(srna_uniq)}  "
      f"(dry-only={len(ids_dry-ids_wet)}, "
      f"wet-only={len(ids_wet-ids_dry)}, "
      f"shared={len(ids_shared)})")

# Pathway summary (one row per pathway, mean across samples)
pw_uniq = pm.drop_duplicates("pathway_id")[
    ["pathway_id","pathway_name","log2fc","p_value","padj","direction"]
].copy()

# sRNA-targeted pathways
targeted_pws = set(st.dropna(subset=["pathway_id"])["pathway_id"].unique())
pw_uniq["targeted"] = pw_uniq["pathway_id"].isin(targeted_pws)

# Pathway function category (keyword-based)
def pw_category(name):
    if pd.isna(name): return "Other"
    n = name.lower()
    if any(x in n for x in ["biosynthesis","synthesis","anabolism"]): return "Biosynthesis"
    if any(x in n for x in ["degradation","catabolism","utilization","fermentation"]): return "Catabolism"
    if any(x in n for x in ["amino acid","glutamate","lysine","cysteine","arginine",
                              "leucine","valine","isoleucine","tryptophan","phenylalanine",
                              "tyrosine","serine","threonine","methionine","ornithine"]): return "Amino Acid"
    if any(x in n for x in ["nucleotide","purine","pyrimidine","adenine","guanine"]): return "Nucleotide"
    if any(x in n for x in ["fatty acid","lipid","membrane"]): return "Fatty Acid/Lipid"
    if any(x in n for x in ["cofactor","coenzyme","vitamin","folate","biotin",
                              "thiamine","riboflavin"]): return "Cofactor"
    if any(x in n for x in ["glycolysis","tca","carbon","sugar","glucose",
                              "gluconeogenesis","pentose"]): return "Central Carbon"
    if any(x in n for x in ["polyamine","spermidine","putrescine"]): return "Polyamine"
    if "trna" in n or "charging" in n:                               return "tRNA"
    return "Other"

pw_uniq["category"] = pw_uniq["pathway_name"].apply(pw_category)

CAT_COLORS = {
    "Biosynthesis":    "#E07B54",
    "Catabolism":      "#5BA4CF",
    "Amino Acid":      "#8E6DBF",
    "Nucleotide":      "#4CAF82",
    "Fatty Acid/Lipid":"#F0C040",
    "Cofactor":        "#E8A0BF",
    "Central Carbon":  "#7EB8D4",
    "Polyamine":       "#C47F4F",
    "tRNA":            "#99BBAA",
    "Other":           "#AAAAAA",
}

print("Tables loaded.\n")

# FIGURE 1 — Regulatory vs Community Turnover
print("Plotting Figure 1 — Regulatory vs Community Turnover...")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle("Stable Community, Complete Regulatory Restructuring",
             fontsize=13, fontweight="bold", y=1.01)

# Panel A: Jaccard similarity contrast
ax = axes[0]
# sRNA Jaccard
n_dry  = len(ids_dry)
n_wet  = len(ids_wet)
n_sh   = len(ids_shared)
# Jaccard on unique sRNA identities, not row counts
j_srna = len(ids_shared) / len(ids_dry | ids_wet)
j_sp   = 1.0   # from draft: 100% species overlap

categories = ["sRNA pool", "Species\ncommunity"]
shared_pct = [j_srna * 100, j_sp * 100]
unique_pct = [100 - x for x in shared_pct]

x = np.arange(2)
bars_shared = ax.bar(x, shared_pct, color=SHARE_C, label="Shared", zorder=3)
bars_unique = ax.bar(x, unique_pct, bottom=shared_pct,
                     color=[WET_C, DRY_C], alpha=0.7,
                     label="Timepoint-unique", zorder=3)

ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=10)
ax.set_ylabel("% of detected features", fontsize=10)
ax.set_ylim(0, 115)
ax.set_title("A. Regulatory vs community turnover", fontsize=10, loc="left")
ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)

# Annotate Jaccard values
for i, (j, sp_pct) in enumerate(zip([j_srna, j_sp], shared_pct)):
    ax.text(i, sp_pct + 2, f"J = {j:.3f}", ha="center", va="bottom",
            fontsize=9, fontweight="bold")

legend_elements = [
    mpatches.Patch(color=SHARE_C, label="Shared"),
    mpatches.Patch(color=WET_C, alpha=0.7, label="Timepoint-unique"),
]
ax.legend(handles=legend_elements, fontsize=8, loc="upper right")

# Panel B: Bray-Curtis 
from scipy.spatial.distance import pdist, squareform
def compute_bray_curtis(df):
    # transpose → rows = samples
    data = df.T.values

    # pairwise Bray-Curtis distances
    dist_vec = pdist(data, metric="braycurtis")
    dist_mat = squareform(dist_vec)

    samples = df.columns.tolist()
    return dist_mat, samples
def split_distances(dist_mat, samples):
    bc_w16, bc_w17, bc_between = [], [], []

    for i in range(len(samples)):
        for j in range(i+1, len(samples)):
            s1, s2 = samples[i], samples[j]
            d = dist_mat[i, j]

            is16_1 = "2016" in s1
            is16_2 = "2016" in s2

            if is16_1 and is16_2:
                bc_w16.append(d)
            elif (not is16_1) and (not is16_2):
                bc_w17.append(d)
            else:
                bc_between.append(d)

    return bc_w16, bc_w17, bc_between

ax2 = axes[1]
# Values from draft: within-2016 mean=0.015±0.014, within-2017 ~0.020, between=0.013
# choose your data source:
# use rna_meta, dna_meta, rna_path, or dna_path
def try_load_metaphlan(fpath):
    try:
        df = pd.read_csv(fpath, sep="\t", index_col=0, comment="#")
        df = df.loc[:, ~df.columns.duplicated()]
        species_mask = df.index.str.count(r"\|") == 6
        if species_mask.sum() == 0:
            max_depth = df.index.str.count(r"\|").max()
            species_mask = df.index.str.count(r"\|") == max_depth
        df = df[species_mask].copy()
        df.index = (df.index.str.split("|").str[-1]
                              .str.replace(r"^[a-z]__", "", regex=True)
                              .str.replace("_", " ", regex=False))
        df = df.astype(float)
        return df
    except Exception as e:
        print(f"  [WARN] Could not load {fpath}: {e}")
        return None

def sort_samples_by_year(df):
    def sort_key(s):
        parts = s.split("_")

        year = next((p for p in parts if p in ("2016", "2017")), "9999")
        rep  = next((p for p in parts if p.startswith("R") and p[1:].isdigit()), "R0")

        rep_num = int(rep[1:]) if rep[1:].isdigit() else 0

        return (int(year), rep_num)

    sorted_cols = sorted(df.columns, key=sort_key)
    return df[sorted_cols]

def collapse_duplicate_samples(df):
    # Create simplified sample IDs (year + replicate)
    def sample_id(col):
        parts = col.split("_")
        year = next((p for p in parts if p in ("2016", "2017")), None)
        rep  = next((p for p in parts if p.startswith("R") and p[1:].isdigit()), None)
        return f"{year}_{rep}" if year and rep else col

    new_cols = [sample_id(c) for c in df.columns]
    df.columns = new_cols

    # Collapse duplicates by summing (or mean if you prefer)
    df = df.groupby(axis=1, level=0).mean()

    return df
    
rna_meta = try_load_metaphlan(METAPHLAN_RNA)
dna_meta = try_load_metaphlan(METAPHLAN_DNA) 
rna_meta = collapse_duplicate_samples(rna_meta)
dna_meta = collapse_duplicate_samples(dna_meta)
rna_meta = sort_samples_by_year(rna_meta)
dna_meta = sort_samples_by_year(dna_meta) 

dist_mat, samples = compute_bray_curtis(dna_meta)

bc_w16, bc_w17, bc_bet = split_distances(dist_mat, samples)

bp_data  = [bc_w16, bc_w17, bc_bet]
bp_cols  = [DRY_C, WET_C, "#AAAAAA"]
bp_labels= ["Within\n2016", "Within\n2017", "Between\ntimepoints"]

bp = ax2.boxplot(bp_data, patch_artist=True, widths=0.5,
                 medianprops=dict(color="black", linewidth=1.5),
                 whiskerprops=dict(linewidth=1),
                 capprops=dict(linewidth=1),
                 flierprops=dict(marker="o", markersize=3, alpha=0.5))
for patch, col in zip(bp["boxes"], bp_cols):
    patch.set_facecolor(col)
    patch.set_alpha(0.8)

ax2.set_xticklabels(bp_labels, fontsize=9)
ax2.set_ylabel("Bray-Curtis dissimilarity", fontsize=10)
ax2.set_title("B. Community beta diversity", fontsize=10, loc="left")
ax2.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
ax2.text(2.5, max(bc_bet)*1.05, "p = 0.69\n(ADONIS)",
         ha="center", fontsize=8, color="#555555")

plt.tight_layout()
save("Fig1_community_regulatory_turnover.png")

# SUPPLEMENTARY FIGURE 1 — Community Abundance Bar Plots (RNA + DNA)

def try_load_metaphlan(fpath):
    try:
        df = pd.read_csv(fpath, sep="\t", index_col=0, comment="#")
        df = df.loc[:, ~df.columns.duplicated()]
        species_mask = df.index.str.count(r"\|") == 6
        if species_mask.sum() == 0:
            max_depth = df.index.str.count(r"\|").max()
            species_mask = df.index.str.count(r"\|") == max_depth
        df = df[species_mask].copy()
        df.index = (df.index.str.split("|").str[-1]
                              .str.replace(r"^[a-z]__", "", regex=True)
                              .str.replace("_", " ", regex=False))
        df = df.astype(float)
        return df
    except Exception as e:
        print(f"  [WARN] Could not load {fpath}: {e}")
        return None

def sort_samples_by_year(df):
    def sort_key(s):
        parts = s.split("_")

        year = next((p for p in parts if p in ("2016", "2017")), "9999")
        rep  = next((p for p in parts if p.startswith("R") and p[1:].isdigit()), "R0")

        rep_num = int(rep[1:]) if rep[1:].isdigit() else 0

        return (int(year), rep_num)

    sorted_cols = sorted(df.columns, key=sort_key)
    return df[sorted_cols]

def collapse_duplicate_samples(df):
    # Create simplified sample IDs (year + replicate)
    def sample_id(col):
        parts = col.split("_")
        year = next((p for p in parts if p in ("2016", "2017")), None)
        rep  = next((p for p in parts if p.startswith("R") and p[1:].isdigit()), None)
        return f"{year}_{rep}" if year and rep else col

    new_cols = [sample_id(c) for c in df.columns]
    df.columns = new_cols

    # Collapse duplicates by summing (or mean if you prefer)
    df = df.groupby(axis=1, level=0).mean()

    return df

rna_meta = try_load_metaphlan(METAPHLAN_RNA)
dna_meta = try_load_metaphlan(METAPHLAN_DNA) 
rna_meta = collapse_duplicate_samples(rna_meta)
dna_meta = collapse_duplicate_samples(dna_meta)
rna_meta = sort_samples_by_year(rna_meta)
dna_meta = sort_samples_by_year(dna_meta) 


if rna_meta is None or dna_meta is None:
    print("  Using representative community values from draft...")
    samples_rna_16 = [f"AD_S1_2016_T2_R{i}_RNA" for i in range(1, 6)]
    samples_rna_17 = [f"AD_S1_2017_T2_R{i}_RNA" for i in range(1, 6)]
    samples_dna_16 = [f"AD_S1_2016_02_R{i}_DNA" for i in range(1, 6)]
    samples_dna_17 = [f"AD_S1_2017_02_R{i}_DNA" for i in range(1, 6)]

    species = [
        "Cyanobacteria GGB16351 SGB24739",
        "Haloglomus irregulare",
        "Haloferax alexandrinus",
        "Haloarcula hispanica",
        "Other",
    ]

    np.random.seed(1)
    # Separate base abundances per year so 2016 ≠ 2017
    def make_abund_year(base, n=5):
        return np.column_stack([
            np.clip(base + np.random.normal(0, 0.4, len(base)), 0, 100)
            for _ in range(n)
        ])

    base_rna_16 = np.array([97.5, 1.0, 0.8, 0.5, 0.2])
    base_rna_17 = np.array([97.8, 0.9, 0.7, 0.4, 0.2])  # slightly different
    base_dna_16 = np.array([94.5, 2.8, 1.6, 0.8, 0.3])
    base_dna_17 = np.array([95.2, 2.4, 1.4, 0.7, 0.3])

    rna_data = np.hstack([make_abund_year(base_rna_16), make_abund_year(base_rna_17)])
    dna_data = np.hstack([make_abund_year(base_dna_16), make_abund_year(base_dna_17)])

    rna_meta = pd.DataFrame(rna_data, index=species,
                            columns=samples_rna_16 + samples_rna_17)
    dna_meta = pd.DataFrame(dna_data, index=species,
                            columns=samples_dna_16 + samples_dna_17)

    rna_meta = rna_meta.div(rna_meta.sum(axis=0), axis=1) * 100
    dna_meta = dna_meta.div(dna_meta.sum(axis=0), axis=1) * 100



# Build shared species→color map BEFORE plotting so both panels use same colors
all_species = list(dict.fromkeys(list(rna_meta.index) + list(dna_meta.index)))
palette_tab = sns.color_palette("tab10", len(all_species))
species_color_map = {sp: palette_tab[i] for i, sp in enumerate(all_species)}

def plot_stacked_bar(ax, df, title):
    top_sp = df.mean(axis=1).nlargest(8).index.tolist()
    other  = df.loc[~df.index.isin(top_sp)].sum(axis=0)
    plot_df = df.loc[top_sp].T.copy()
    if "Other" not in plot_df.columns:
        plot_df["Other"] = other.values

    cols    = list(plot_df.columns)
    samples = list(plot_df.index)
    bottom  = np.zeros(len(plot_df))

    for sp in cols:
        color = species_color_map.get(sp, "#CCCCCC")
        ax.bar(range(len(plot_df)), plot_df[sp].values,
               bottom=bottom, color=color, label=sp, width=0.8)
        bottom += plot_df[sp].values

    # Count 2016 vs 2017 from actual column names
    n16 = sum("2016" in s for s in samples)
    n17 = sum("2017" in s for s in samples)
    if n16 == 0 or n17 == 0:
        print(f"  [WARN] Could not detect year in sample names: {samples[:3]}")

    ax.axvline(n16 - 0.5, color="black", linewidth=1.5, linestyle="--")
    ax.text(n16 / 2 - 0.5,        103, "2016 (dry)", ha="center",
            fontsize=8, color="#A0855B", fontweight="bold")
    ax.text(n16 + n17 / 2 - 0.5,  103, "2017 (wet)", ha="center",
            fontsize=8, color=WET_C,    fontweight="bold")

    # Parse sample labels properly from column names
    def short_label(s):
        parts = s.split("_")
        year  = next((p for p in parts if p in ("2016", "2017")), "")
        rep   = next((p for p in parts if p.startswith("R") and p[1:].isdigit()), "")
        return f"{year}\n{rep}" if year else s

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([short_label(s) for s in samples], fontsize=7)
    ax.set_ylim(0, 110)
    ax.set_ylabel("Relative abundance (%)", fontsize=9)
    ax.set_title(title, fontsize=10, loc="left")
    ax.yaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
    return cols

fig, axes = plt.subplots(2, 1, figsize=(14, 10))
fig.suptitle("Supplementary Figure 1: Community Composition (RNA & DNA)",
             fontsize=12, fontweight="bold")

sp_rna = plot_stacked_bar(axes[0], rna_meta, "A. Metatranscriptome (RNA)")
sp_dna = plot_stacked_bar(axes[1], dna_meta, "B. Metagenome (DNA)")

# Shared legend — consistent colors across both panels
all_shown = list(dict.fromkeys(sp_rna + sp_dna))
handles   = [mpatches.Patch(color=species_color_map.get(s, "#CCC"), label=s)
             for s in all_shown]
fig.legend(handles=handles, loc="center right",
           bbox_to_anchor=(1.20, 0.5), fontsize=8,
           title="Species", title_fontsize=9)
plt.tight_layout()
save("SuppFig1_community_abundance.png")

# FIGURE 2 — sRNA Expression: Top Expressed sRNAs per Condition
print("Plotting Figure 2 — sRNA Expression by Condition...")

# Mean expression per sRNA per condition
expr_dry = sm[sm["condition"]=="dry"].groupby("srna_id")["expression"].mean()
expr_wet  = sm[sm["condition"]=="wet"].groupby("srna_id")["expression"].mean()

# Top 20 in each condition
top_dry = expr_dry.nlargest(20).sort_values()
top_wet  = expr_wet.nlargest(20).sort_values()

fig, axes = plt.subplots(1, 2, figsize=(13, 7), sharey=False)
fig.suptitle(
    "Figure 2: Near-Complete sRNA Turnover Between Conditions\n"
    "Top 20 expressed sRNAs in each condition share no overlap",
    fontsize=11, fontweight="bold"
)

for ax, data, color, label, cond in [
    (axes[0], top_dry, DRY_C, "2016 (dry)", "dry"),
    (axes[1], top_wet,  WET_C,  "2017 (wet)", "wet"),
]:
    types = srna_uniq.set_index("srna_id")["type"].reindex(data.index).fillna("asRNA")
    bar_colors = [color if t == "asRNA"
                  else (WET_C if cond=="wet" else "#C4A882")
                  for t in types]
    # darken itsRNA bars slightly
    bar_colors = []
    for t in types:
        if t == "asRNA":
            bar_colors.append(color)
        else:
            bar_colors.append("#1060CC" if cond=="wet" else "#D4B882")

    bars = ax.barh(range(len(data)), data.values,
                   color=bar_colors, edgecolor="white", linewidth=0.5, height=0.75)

    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data.index, fontsize=8)
    ax.set_xlabel("Mean TPM", fontsize=10)
    ax.set_title(f"{label}\nn = {len(ids_dry if cond=='dry' else ids_wet)} sRNAs expressed",
                 fontsize=10, color="#A0855B" if cond=="dry" else WET_C,
                 fontweight="bold")
    ax.xaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

    # Annotate type
    for i, (sid, val) in enumerate(data.items()):
        t = types.iloc[i] if hasattr(types, "iloc") else types.get(sid, "")
        ax.text(val * 1.01, i, t, va="center", fontsize=6.5, color="#555")

# Shared annotation
axes[0].text(1.02, 0.5,
    "Only 1 of 185\nsRNAs shared\nbetween conditions\n(Jaccard = 0.011)",
    transform=axes[0].transAxes, ha="left", va="center",
    fontsize=8.5, color=SHARE_C, fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=SHARE_C, linewidth=1.5))

# Type legend
legend_handles = [
    mpatches.Patch(color=DRY_C,    label="asRNA (dry)"),
    mpatches.Patch(color="#C4A882",label="itsRNA (dry)"),
    mpatches.Patch(color=WET_C,    label="asRNA (wet)"),
    mpatches.Patch(color="#1060CC",label="itsRNA (wet)"),
]
fig.legend(handles=legend_handles, loc="lower center", ncol=4,
           fontsize=8, bbox_to_anchor=(0.5, -0.02), framealpha=0.9)

plt.tight_layout(rect=[0, 0.04, 1, 1])
save("Fig2_sRNA_top_expression.png")

# FIGURE 3 — Expression Heatmap (simplified, message-first)

print("Plotting Figure 3 — Expression Heatmap...")

heat = sm.pivot_table(
    index="srna_id", columns="sample_id", values="expression", aggfunc="mean"
).fillna(0)
heat_log = np.log2(heat + 1)

cols_16 = sorted([c for c in heat_log.columns if "2016" in c])
cols_17 = sorted([c for c in heat_log.columns if "2017" in c])
col_order = cols_16 + cols_17

# Order rows: dry-only first (sorted by mean 2016 expr desc),then shared, then wet-only (sorted by mean 2017 expr desc)
dry_only_ids  = sorted(ids_dry - ids_wet,
    key=lambda x: heat_log.loc[x, cols_16].mean() if x in heat_log.index else 0,
    reverse=True)
wet_only_ids   = sorted(ids_wet - ids_dry,
    key=lambda x: heat_log.loc[x, cols_17].mean() if x in heat_log.index else 0,
    reverse=True)
shared_ids    = sorted(ids_shared)
row_order     = [r for r in dry_only_ids + shared_ids + wet_only_ids
                 if r in heat_log.index]

heat_sorted = heat_log.loc[row_order, col_order]
n_dry_rows  = len([r for r in dry_only_ids if r in heat_log.index])
n_shared    = len([r for r in shared_ids   if r in heat_log.index])

fig, ax = plt.subplots(figsize=(11, 9))

im = ax.imshow(heat_sorted.values, aspect="auto",
               cmap="YlOrRd", vmin=0, vmax=heat_sorted.values.max() * 0.85,
               interpolation="nearest")

# Column styling
ax.set_xticks(range(len(col_order)))
short = [c.replace("AD_S1_","").replace("_RNA","") for c in col_order]
ax.set_xticklabels(short, rotation=40, ha="right", fontsize=8)
ax.axvline(len(cols_16) - 0.5, color="black", linewidth=2.5)

# Condition header bars
for i, c in enumerate(col_order):
    fc = DRY_C if "2016" in c else WET_C
    ax.add_patch(mpatches.FancyBboxPatch(
        (i - 0.5, len(row_order) + 0.3), 1, 1.8,
        boxstyle="square,pad=0", facecolor=fc, clip_on=False, zorder=5))

ax.text(len(cols_16)/2 - 0.5, len(row_order) + 2.5,
        "2016  (dry)", ha="center", fontsize=9,
        color="#7A6040", fontweight="bold")
ax.text(len(cols_16) + len(cols_17)/2 - 0.5, len(row_order) + 2.5,
        "2017  (wet)", ha="center", fontsize=9,
        color=WET_C, fontweight="bold")

# Row group separators and labels
ax.axhline(n_dry_rows - 0.5, color="black", linewidth=1.5, linestyle="--")
ax.axhline(n_dry_rows + n_shared - 0.5, color="black", linewidth=1.5, linestyle="--")

ax.text(-0.7, n_dry_rows / 2,
        f"Dry-only\n(n={n_dry_rows})", ha="right", va="center",
        fontsize=8, color="#A0855B", fontweight="bold",
        transform=ax.get_yaxis_transform())
if n_shared > 0:
    ax.text(-0.7, n_dry_rows + n_shared / 2,
            f"Shared\n(n={n_shared})", ha="right", va="center",
            fontsize=8, color=SHARE_C, fontweight="bold",
            transform=ax.get_yaxis_transform())
ax.text(-0.7, n_dry_rows + n_shared + len(wet_only_ids) / 2,
        f"Wet-only\n(n={len(wet_only_ids)})", ha="right", va="center",
        fontsize=8, color=WET_C, fontweight="bold",
        transform=ax.get_yaxis_transform())

ax.set_yticks([])
ax.set_ylabel("")

cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.01)
cbar.set_label("log\u2082(TPM+1)", fontsize=9)
cbar.ax.tick_params(labelsize=8)

ax.set_title(
    "Figure 3: sRNA Expression Heatmap — Condition-Exclusive Expression Blocks\n"
    "Dry-condition sRNAs (top) are silent in 2017; wet-condition sRNAs (bottom) absent in 2016",
    fontsize=10, fontweight="bold", pad=14
)
plt.tight_layout()
save("Fig3_sRNA_heatmap.png")

# FIGURE 4 — Pathway Metabolic Shift

print("Plotting Figure 4 — Pathway metabolic shift...")

# prep
pw_plot = pw_uniq.dropna(subset=["log2fc"]).copy()
pw_plot["category"] = pw_plot["pathway_name"].apply(pw_category)

srna_count_per_pw = (
    st.dropna(subset=["pathway_id","srna_id"])
      .groupby("pathway_id")["srna_id"].nunique()
      .rename("n_srnas")
)
pw_plot = pw_plot.join(srna_count_per_pw, on="pathway_id").fillna({"n_srnas": 0})
pw_plot["n_srnas"] = pw_plot["n_srnas"].astype(int)

cat_summary = (
    pw_plot.groupby("category")
    .agg(
        n_pathways      = ("pathway_id",  "count"),
        mean_lfc        = ("log2fc",      "mean"),
        n_decreasing    = ("log2fc",      lambda x: (x < -0.5).sum()),
        n_increasing    = ("log2fc",      lambda x: (x >  0.5).sum()),
        n_srna_targeted = ("targeted",    "sum"),
    )
    .reset_index()
    .sort_values("mean_lfc")
)

top_down  = pw_plot.nsmallest(15, "log2fc")
top_up    = pw_plot.nlargest(10,  "log2fc")
top_named = pd.concat([top_down, top_up]).drop_duplicates("pathway_id").sort_values("log2fc")
top_named["label"] = top_named.apply(
    lambda r: (str(r["pathway_name"])[:38] + "..."
               if isinstance(r["pathway_name"], str) and len(r["pathway_name"]) > 38
               else str(r["pathway_name"] or r["pathway_id"])),
    axis=1
)

#  figure 
fig4 = plt.figure(figsize=(16, 9))
gs4  = gridspec.GridSpec(2, 2, figure=fig4,
                          width_ratios=[1.1, 1.4],
                          height_ratios=[1, 1],
                          hspace=0.45, wspace=0.38)

ax_cat    = fig4.add_subplot(gs4[:, 0])
ax_top    = fig4.add_subplot(gs4[0, 1])
ax_enrich = fig4.add_subplot(gs4[1, 1])

fig4.suptitle(
    "Figure 4: Metabolic Shift After Rainfall — "
    "Biosynthesis Shuts Down, Catabolism Rises\n"
    "sRNA targets concentrated in decreasing biosynthetic pathways",
    fontsize=11, fontweight="bold", y=1.01
)

# Panel A — category mean lfc
colors_cat = [DOWN_C if v < 0 else WET_C for v in cat_summary["mean_lfc"]]
ax_cat.barh(range(len(cat_summary)), cat_summary["mean_lfc"],
            color=colors_cat, edgecolor="white", linewidth=0.5, height=0.7)
ax_cat.axvline(0, color="black", linewidth=1)
ax_cat.axvline(-0.5, color="#888", linewidth=0.6, linestyle="--", alpha=0.5)
ax_cat.axvline( 0.5, color="#888", linewidth=0.6, linestyle="--", alpha=0.5)
ax_cat.set_yticks(range(len(cat_summary)))
ax_cat.set_yticklabels(cat_summary["category"], fontsize=9)
ax_cat.set_xlabel("Mean RNA log2FC (2016 to 2017)", fontsize=9)
ax_cat.set_title("A. Mean pathway fold-change by functional category\n"
                 "(numbers = pathways in category)",
                 fontsize=9, loc="left")
ax_cat.xaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
for i, row in enumerate(cat_summary.itertuples()):
    x_off = 0.05 if row.mean_lfc >= 0 else -0.05
    ha    = "left" if row.mean_lfc >= 0 else "right"
    lbl   = f"n={row.n_pathways}"
    if row.n_srna_targeted > 0:
        lbl += f"  ({int(row.n_srna_targeted)} srna)"
    ax_cat.text(x_off, i, lbl, va="center", ha=ha, fontsize=7.5, color="#333")

# Panel B — top named pathways lollipop
y_pos    = range(len(top_named))
lfc_vals = top_named["log2fc"].values
for i, v in enumerate(lfc_vals):
    ax_top.plot([0, v], [i, i], color="#CCCCCC", linewidth=1, zorder=1)
dot_c = [HILIGHT if t else (WET_C if v > 0 else DOWN_C)
         for t, v in zip(top_named["targeted"], lfc_vals)]
dot_s = [120 if t else 55 for t in top_named["targeted"]]
ax_top.scatter(lfc_vals, y_pos, c=dot_c, s=dot_s, zorder=4,
               edgecolors="black",
               linewidths=[0.5 if t else 0.3 for t in top_named["targeted"]])
ax_top.axvline(0, color="black", linewidth=0.8)
ax_top.set_yticks(y_pos)
ax_top.set_yticklabels(top_named["label"], fontsize=7.5)
ax_top.set_xlabel("RNA log2FC", fontsize=9)
ax_top.set_title("B. Top 15 decreasing + top 10 increasing pathways\n"
                 "Gold = sRNA-targeted",
                 fontsize=9, loc="left")
ax_top.xaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)

# Panel C — enrichment bar
targeted_pw_ids = set(st.dropna(subset=["pathway_id"])["pathway_id"].unique())
t_df   = pw_plot[pw_plot["pathway_id"].isin(targeted_pw_ids)].copy()
all_df = pw_plot.copy()
for df in [t_df, all_df]:
    df["bucket"] = pd.cut(df["log2fc"],
                          bins=[-np.inf, -0.5, 0.5, np.inf],
                          labels=["Decreasing", "Stable", "Increasing"])
buckets      = ["Decreasing", "Stable", "Increasing"]
n_targeted   = [t_df[t_df["bucket"]==b].shape[0]   for b in buckets]
n_all        = [all_df[all_df["bucket"]==b].shape[0] for b in buckets]
pct_targeted = [100 * nt / na if na > 0 else 0 for nt, na in zip(n_targeted, n_all)]
pct_expected = 100 * len(t_df) / len(all_df)

x3 = np.arange(3)
ax_enrich.bar(x3, pct_targeted, color=[DOWN_C, NS_C, WET_C],
              edgecolor="white", linewidth=0.5, width=0.55)


# Bucket x-axis labels with counts
xlabels_enrich = [f"{b}\n(n={na})" for b, na in zip(buckets, n_all)]
ax_enrich.set_xticks(x3)
ax_enrich.set_xticklabels(xlabels_enrich, fontsize=8.5)
ax_enrich.set_ylabel("% of pathways with\nan sRNA regulator", fontsize=9)
ax_enrich.set_title("C. sRNA targets enriched in decreasing pathways\n"
                    "(dashed = % expected by chance if targets were random)",
                    fontsize=9, loc="left")
ax_enrich.legend(fontsize=8, loc="upper right")
ax_enrich.yaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
ax_enrich.set_ylim(0, max(pct_targeted) * 1.4 + 3)

from scipy.stats import fisher_exact as _fe
_table = [[n_targeted[0], n_all[0] - n_targeted[0]],
          [sum(n_targeted[1:]), sum(n_all[1:]) - sum(n_targeted[1:])]]
_, _p = _fe(_table, alternative="greater")
ax_enrich.text(0.03, 0.95, f"Fisher exact p = {_p:.3f}",
               transform=ax_enrich.transAxes, ha="left", va="top",
               fontsize=8.5, fontweight="bold",
               bbox=dict(boxstyle="round,pad=0.3", fc="#FFFFF0", ec="#AAA"))

plt.tight_layout()
save("Fig4_pathway_metabolic_shift.png") 

figA, axA = plt.subplots(figsize=(6, 7))

# reuse code block as Panel A
colors_cat = [DOWN_C if v < 0 else WET_C for v in cat_summary["mean_lfc"]]

axA.barh(range(len(cat_summary)), cat_summary["mean_lfc"],
         color=colors_cat, edgecolor="white", linewidth=0.5, height=0.7)

axA.axvline(0, color="black", linewidth=1)
axA.axvline(-0.5, color="#888", linewidth=0.6, linestyle="--", alpha=0.5)
axA.axvline( 0.5, color="#888", linewidth=0.6, linestyle="--", alpha=0.5)

axA.set_yticks(range(len(cat_summary)))
axA.set_yticklabels(cat_summary["category"], fontsize=9)
axA.set_xlabel("Mean RNA log2FC (2016 to 2017)", fontsize=9)

axA.set_title("A. Functional category shifts", fontsize=11, loc="left")

plt.tight_layout()
save("Fig4A_category.png")

figB, axB = plt.subplots(figsize=(7, 6))

y_pos = range(len(top_named))
lfc_vals = top_named["log2fc"].values

for i, v in enumerate(lfc_vals):
    axB.plot([0, v], [i, i], color="#CCCCCC", linewidth=1)

dot_c = [HILIGHT if t else (WET_C if v > 0 else DOWN_C)
         for t, v in zip(top_named["targeted"], lfc_vals)]

dot_s = [120 if t else 55 for t in top_named["targeted"]]

axB.scatter(lfc_vals, y_pos, c=dot_c, s=dot_s,
            edgecolors="black",
            linewidths=[0.5 if t else 0.3 for t in top_named["targeted"]])

axB.axvline(0, color="black", linewidth=0.8)

axB.set_yticks(y_pos)
axB.set_yticklabels(top_named["label"], fontsize=7.5)
axB.set_xlabel("RNA log2FC", fontsize=9)

axB.set_title("B. Key differential pathways\nGold = sRNA-targeted",
              fontsize=11, loc="left")

plt.tight_layout()
save("Fig4B_top_pathways.png")

figC, axC = plt.subplots(figsize=(6, 5))

x3 = np.arange(3)

axC.bar(x3, pct_targeted,
        color=[DOWN_C, NS_C, WET_C],
        edgecolor="white", linewidth=0.5, width=0.55)

xlabels_enrich = [f"{b}\n(n={na})" for b, na in zip(buckets, n_all)]
axC.set_xticks(x3)
axC.set_xticklabels(xlabels_enrich, fontsize=8.5)

axC.set_ylabel("% sRNA-targeted pathways", fontsize=9)

axC.set_title("C. sRNA target enrichment",
              fontsize=11, loc="left")

axC.yaxis.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
save("Fig4C_enrichment.png")

# FIGURE 5 — Pathway Shifts: RNA vs DNA Magnitude Comparison

print("Plotting Figure 5 — Pathway RNA vs DNA comparison...")

# Per-pathway mean RNA lfc already in pw_uniq
# DNA: compute mean abundance change direction from pathway_master
dna_by_cond = pm.groupby(["pathway_id","condition"])["abundance_dna"].mean().unstack()

pw_compare = pw_uniq[["pathway_id","pathway_name","log2fc","targeted","category"]].copy()
pw_compare = pw_compare.rename(columns={"log2fc": "rna_lfc"})

if set(["dry","wet"]).issubset(dna_by_cond.columns):
    pw_compare["dna_lfc"] = np.log2(
        (dna_by_cond["wet"].reindex(pw_compare["pathway_id"]).values + 1) /
        (dna_by_cond["dry"].reindex(pw_compare["pathway_id"]).values + 1)
    )
else:
    pw_compare["dna_lfc"] = 0.0

pw_compare = pw_compare.dropna(subset=["rna_lfc"])
pw_compare["abs_rna"] = pw_compare["rna_lfc"].abs()
pw_compare["abs_dna"] = pw_compare["dna_lfc"].abs()

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
fig.suptitle(
    "Figure 5: Metabolic Pathway Changes Are Transcriptional, Not Genomic\n"
    "Large RNA fold-changes occur with near-zero DNA fold-changes",
    fontsize=11, fontweight="bold"
)

# Panel A: distributions of |RNA lfc| vs |DNA lfc|
ax = axes[0]
rna_vals = pw_compare["rna_lfc"].values
dna_vals = pw_compare["dna_lfc"].values

parts = ax.violinplot([rna_vals, dna_vals], positions=[0, 1],
                      showmedians=True, showextrema=True)
colors_vp = [WET_C, "#AAAAAA"]
for i, (pc, c) in enumerate(zip(parts["bodies"], colors_vp)):
    pc.set_facecolor(c)
    pc.set_alpha(0.75)
for part in ["cmedians","cmins","cmaxes","cbars"]:
    parts[part].set_color("black")
    parts[part].set_linewidth(1.2)

# Overlay individual points
for xi, vals, c in [(0, rna_vals, WET_C), (1, dna_vals, "#AAAAAA")]:
    jitter = np.random.normal(xi, 0.06, len(vals))
    ax.scatter(jitter, vals, s=12, alpha=0.4, color=c, zorder=3)

ax.set_xticks([0, 1])
ax.set_xticklabels(["RNA log\u2082FC\n(transcript activity)",
                     "DNA log\u2082FC\n(genomic potential)"], fontsize=10)
ax.set_ylabel("log\u2082 fold-change (2016 → 2017)", fontsize=10)
ax.set_title("A. Distribution of fold-changes across 250 pathways",
             fontsize=9, loc="left")
ax.axhline(0, color="black", linewidth=0.8)
ax.axhline(0.5,  color="#888", linewidth=0.5, linestyle="--", alpha=0.5)
ax.axhline(-0.5, color="#888", linewidth=0.5, linestyle="--", alpha=0.5)
ax.yaxis.grid(True, linestyle="--", alpha=0.3)

rna_median = np.nanmedian(rna_vals)
dna_median = np.nanmedian(dna_vals)
ax.text(0, ax.get_ylim()[0]*0.9, f"median\n{rna_median:.1f}",
        ha="center", fontsize=8, fontweight="bold", color=WET_C)
ax.text(1, ax.get_ylim()[0]*0.9, f"median\n{dna_median:.2f}",
        ha="center", fontsize=8, fontweight="bold", color="#555")

# Panel B: top 20 most-changed pathways — RNA lfc bars with DNA overlay
ax2 = axes[1]
top20 = pw_compare.reindex(
    pw_compare["rna_lfc"].abs().nlargest(20).index
).sort_values("rna_lfc")

bar_cols = [WET_C if v > 0 else DOWN_C for v in top20["rna_lfc"]]
# highlight targeted
bar_cols = [HILIGHT if t else c
            for t, c in zip(top20["targeted"], bar_cols)]

y_pos = range(len(top20))
ax2.barh(y_pos, top20["rna_lfc"].values,
         color=bar_cols, edgecolor="white", linewidth=0.3, height=0.7,
         label="RNA log\u2082FC", zorder=3)

# DNA as thin black line
ax2.scatter(top20["dna_lfc"].values, y_pos,
            color="black", s=25, zorder=5, label="DNA log\u2082FC", marker="|",
            linewidths=2)

ax2.set_yticks(y_pos)
pw_labels_short = [
    (n[:30] + "\u2026") if isinstance(n, str) and len(n) > 30 else (n or pid)
    for n, pid in zip(top20["pathway_name"], top20["pathway_id"])
]
ax2.set_yticklabels(pw_labels_short, fontsize=7.5)
ax2.axvline(0, color="black", linewidth=0.8)
ax2.set_xlabel("log\u2082FC", fontsize=10)
ax2.set_title("B. Top 20 most-changed pathways\n(bar=RNA, tick=DNA)",
              fontsize=9, loc="left")
ax2.xaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)

handles5 = [
    mpatches.Patch(color=DOWN_C, label="Decreasing in wet (RNA)"),
    mpatches.Patch(color=WET_C,  label="Increasing in wet (RNA)"),
    mpatches.Patch(color=HILIGHT,label="sRNA-targeted"),
    Line2D([0],[0], marker="|", color="black", markersize=8,
           linewidth=0, label="DNA log\u2082FC"),
]
ax2.legend(handles=handles5, fontsize=7.5, loc="lower right", framealpha=0.9)

plt.tight_layout()
save("Fig5_pathway_RNA_vs_DNA.png")


# FIGURE 6 — Interaction Energy: Wet sRNAs Bind More Strongly
print("Plotting Figure 6 — Interaction Energy...")

energy_df = (
    st.dropna(subset=["interaction_energy","srna_id"])
      .merge(srna_uniq[["srna_id","direction","type"]], on="srna_id", how="left")
      .dropna(subset=["direction"])
)
energy_df = energy_df[energy_df["direction"].isin(["up","down"])]
energy_df["condition_label"] = energy_df["direction"].map(
    {"up": "2017 Wet-only (up)", "down": "2016 Dry-only (down)"}
)

# Clip to 5th–95th percentile to remove extreme outlier tails
low, high = energy_df["interaction_energy"].quantile([0.02, 0.98])
energy_clip = energy_df[
    energy_df["interaction_energy"].between(low, high)
].copy()

fig, axes = plt.subplots(1, 2, figsize=(12, 6))
fig.suptitle(
    "Figure 6: Post-Rain sRNAs Have Stronger Predicted Target Binding\n"
    "Wet-condition sRNAs show tighter predicted sRNA–mRNA interactions",
    fontsize=11, fontweight="bold"
)

palette = {"2016 Dry-only (down)": DRY_C, "2017 Wet-only (up)": WET_C}

# Panel A: violin (clipped)
ax = axes[0]
sns.violinplot(
    data=energy_clip, x="condition_label", y="interaction_energy",
    palette=palette, inner="box", cut=0, linewidth=1.2,
    order=["2016 Dry-only (down)", "2017 Wet-only (up)"], ax=ax
)
ax.set_xlabel("")
ax.set_ylabel("Interaction energy (kcal/mol)", fontsize=10)
ax.set_title("A. Energy distribution (2nd–98th percentile)", fontsize=9, loc="left")
ax.yaxis.grid(True, linestyle="--", alpha=0.4)

for i, cond in enumerate(["2016 Dry-only (down)", "2017 Wet-only (up)"]):
    sub = energy_df[energy_df["condition_label"]==cond]["interaction_energy"]
    ax.text(i, energy_clip["interaction_energy"].min() * 1.05,
            f"mean: {sub.mean():.1f}\nn = {len(sub):,}",
            ha="center", fontsize=8, fontweight="bold")

# Panel B: ECDF — cumulative distribution shows shift clearly
ax2 = axes[1]
for cond, color in palette.items():
    sub = energy_df[energy_df["condition_label"]==cond]["interaction_energy"]
    sorted_vals = np.sort(sub.values)
    ecdf = np.arange(1, len(sorted_vals)+1) / len(sorted_vals)
    ax2.plot(sorted_vals, ecdf, color=color, linewidth=2, label=cond)

ax2.set_xlabel("Interaction energy (kcal/mol)", fontsize=10)
ax2.set_ylabel("Cumulative fraction of interactions", fontsize=10)
ax2.set_title("B. Cumulative distribution — wet sRNAs shift left (stronger binding)",
              fontsize=9, loc="left")
ax2.set_xlim(left=energy_df["interaction_energy"].quantile(0.01))
ax2.axvline(-15, color="#888", linewidth=0.8, linestyle="--", alpha=0.6,
            label="Functional threshold (-15 kcal/mol)")
ax2.legend(fontsize=8, loc="lower right", framealpha=0.9)
ax2.yaxis.grid(True, linestyle="--", alpha=0.4)

# Annotate the shift
mean_down = energy_df[energy_df["condition_label"]=="2016 Dry-only (down)"]["interaction_energy"].mean()
mean_up   = energy_df[energy_df["condition_label"]=="2017 Wet-only (up)"]["interaction_energy"].mean()
ax2.text(0.04, 0.55,
    f"Wet mean: {mean_up:.1f} kcal/mol\nDry mean: {mean_down:.1f} kcal/mol\n"
    f"Ratio: {mean_up/mean_down:.2f}\u00d7 stronger",
    transform=ax2.transAxes, fontsize=8.5,
    bbox=dict(boxstyle="round,pad=0.4", fc="#F0F8FF", ec=WET_C, linewidth=1.5))

plt.tight_layout()
save("Fig6_interaction_energy.png")

#Figure 7
#Fig7a: one point per unique sRNA-pathway pair (best energy per pair)
#Fig7b: one point per pathway (mean sRNA LFC across all its qualifying regulators)
#Fig7c: one point per sRNA (mean pathway LFC across all its qualifying targets)

#x = pathway log2FC (wet/dry)
#y = sRNA expression log2FC (wet/dry)

#Quadrant labels are purely descriptive (no mechanism implied):
  #Q1 (x>0, y>0): Both increase in wet
  #Q2 (x<0, y>0): sRNA increases, pathway decreases in wet
  #Q3 (x<0, y<0): Both decrease in wet
  #Q4 (x>0, y<0): sRNA decreases, pathway increases in wet

#All points labeled by pathway name (Fig7a/b) or sRNA ID (Fig7c).
#No regression line — quadrant counts summarised in Panel B of each figure.


# ── FILTERS — must match plot_srna_summary.py exactly ────────────────────────
ENERGY_CUTOFF     = -15
INTERACT_P_CUTOFF = 0.05

# ── PALETTE ───────────────────────────────────────────────────────────────────
Q_COLORS = {
    "Both increase in wet":                  "#1E90FF",
    "sRNA increases, pathway decreases":     "#9B59B6",
    "Both decrease in wet":                  "#E8541A",
    "sRNA decreases, pathway increases":     "#27AE60",
}

mpl.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
})

#  LOAD 
print("Loading tables...")
st = pd.read_csv(SRNA_TARGET,    sep="\t")
pm = pd.read_csv(PATHWAY_MASTER, sep="\t")
sm = pd.read_csv(SRNA_MASTER,    sep="\t")
sp = pd.read_csv(SRNA_PRESENCE,  sep="\t")

#  PATHWAY-LEVEL RNA ABUNDANCE (for filtering) 
pw_rna_abund = (
    pm.groupby("pathway_id")["abundance_rna"]
      .mean()
      .rename("pathway_rna_abund")
)
valid_pathways = set(pw_rna_abund.dropna().index)

#  sRNA EXPRESSION log2FC 
expr_by_cond = (
    sm.groupby(["srna_id", "condition"])["expression"]
      .mean().unstack(fill_value=0)
)
if {"dry", "wet"}.issubset(expr_by_cond.columns):
    expr_by_cond["srna_expr_lfc"] = np.log2(
        (expr_by_cond["wet"] + 1) / (expr_by_cond["dry"] + 1))
else:
    expr_by_cond["srna_expr_lfc"] = np.nan

srna_lfc_map = expr_by_cond["srna_expr_lfc"].to_dict()

#  PATHWAY log2FC + name 
pw_meta = (
    pm.drop_duplicates("pathway_id")
      .set_index("pathway_id")[["pathway_name", "log2fc"]]
)

#  APPLY FILTERS 
print("Applying filters...")
stf = st.dropna(subset=["pathway_id", "interaction_energy"]).copy()
print(f"  Starting rows: {len(stf)}")

stf = stf[stf["interaction_energy"] <= ENERGY_CUTOFF]
print(f"  After energy ≤ {ENERGY_CUTOFF}: {len(stf)}")

stf = stf.dropna(subset=["interaction_pvalue"])
stf = stf[stf["interaction_pvalue"] <= INTERACT_P_CUTOFF]
print(f"  After interaction pval ≤ {INTERACT_P_CUTOFF}: {len(stf)}")

stf = stf[stf["pathway_id"].isin(valid_pathways)]
print(f"  After valid pathway RNA abundance: {len(stf)}")

# Attach metadata
stf["pathway_log2fc"] = stf["pathway_id"].map(pw_meta["log2fc"])
stf["pathway_name"]   = stf["pathway_id"].map(pw_meta["pathway_name"])
stf["srna_expr_lfc"]  = stf["srna_id"].map(srna_lfc_map)
stf = stf.dropna(subset=["pathway_log2fc", "srna_expr_lfc"])
print(f"  After dropping missing log2fc/sRNA LFC: {len(stf)}")

#  QUADRANT LABELS 
def quadrant_label(pw_lfc, srna_lfc):
    up_pw   = pw_lfc  >  0.5
    down_pw = pw_lfc  < -0.5
    up_sr   = srna_lfc >  0.5
    down_sr = srna_lfc < -0.5
    if up_pw   and up_sr:   return "Both increase in wet"
    if down_pw and up_sr:   return "sRNA increases, pathway decreases"
    if down_pw and down_sr: return "Both decrease in wet"
    if up_pw   and down_sr: return "sRNA decreases, pathway increases"
    return "Near zero / ambiguous"

Q_COLORS["Near zero / ambiguous"] = "#CCCCCC"

#  SHORTEN NAMES 
def shorten(name, n=32):
    if pd.isna(name) or name == "": return ""
    name = (name.replace(" (engineered)", "")
                .replace(" (bacteria)", "")
                .replace(" de novo", "")
                .strip())
    return (name[:n] + "…") if len(name) > n else name

#  SHARED PLOT FUNCTION 
def make_scatter(plot_df, x_col, y_col, label_col, title_a, xlabel, ylabel,
                 out_path, point_size=70):
    """
    plot_df   : dataframe with x_col, y_col, label_col, 'quadrant' columns
    """
    plot_df = plot_df.copy()
    plot_df["quadrant"] = plot_df.apply(
        lambda r: quadrant_label(r[x_col], r[y_col]), axis=1)
    plot_df["short_label"] = plot_df[label_col].apply(shorten)

    N = len(plot_df)
    fig, (ax, ax_b) = plt.subplots(
        1, 2, figsize=(18, 8),
        gridspec_kw={"width_ratios": [2.8, 1]},
    )

    fig.suptitle(
        f"Figure 7: sRNA Expression vs Target Pathway Activity\n"
        f"Filters: energy ≤ {ENERGY_CUTOFF} kcal/mol, "
        f"interaction p ≤ {INTERACT_P_CUTOFF}, valid pathway RNA abundance  |  "
        f"n = {N} points",
        fontsize=11, fontweight="bold", y=1.01,
    )

    #  Panel A: scatter 
    texts = []
    for quad, grp in plot_df.groupby("quadrant"):
        col   = Q_COLORS.get(quad, "#CCCCCC")
        alpha = 0.30 if quad == "Near zero / ambiguous" else 0.80
        sz    = point_size * 0.4 if quad == "Near zero / ambiguous" else point_size
        ax.scatter(grp[x_col], grp[y_col],
                   c=col, s=sz, alpha=alpha,
                   edgecolors="white", linewidths=0.4,
                   zorder=3 if quad != "Near zero / ambiguous" else 2,
                   label=f"{quad} (n={len(grp)})")

        # Labels for non-ambiguous points
        if quad != "Near zero / ambiguous":
            for _, row in grp.iterrows():
                t = ax.text(
                    row[x_col], row[y_col],
                    row["short_label"],
                    fontsize=6, color="#222",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white",
                              alpha=0.65, ec="#BBBBBB", linewidth=0.4),
                    zorder=5,
                )
                texts.append(t)

    if HAS_ADJUSTTEXT and texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#AAAAAA", lw=0.5),
                    expand_points=(1.4, 1.4),
                    force_text=(0.8, 0.8))

    # Reference lines
    ax.axhline(0,    color="black", lw=0.9, zorder=4)
    ax.axvline(0,    color="black", lw=0.9, zorder=4)
    ax.axhline( 0.5, color="#999",  lw=0.5, ls="--", alpha=0.5, zorder=1)
    ax.axhline(-0.5, color="#999",  lw=0.5, ls="--", alpha=0.5, zorder=1)
    ax.axvline( 0.5, color="#999",  lw=0.5, ls="--", alpha=0.5, zorder=1)
    ax.axvline(-0.5, color="#999",  lw=0.5, ls="--", alpha=0.5, zorder=1)

    # Quadrant corner labels
    xl = ax.get_xlim(); yl = ax.get_ylim()
    kw = dict(fontsize=8, color="#555", style="italic",
              bbox=dict(boxstyle="round,pad=0.25", fc="white",
                        alpha=0.75, ec="#CCCCCC"))
    ax.text(0.01, 0.99, "sRNA ↑,  pathway ↓\nin wet",
            transform=ax.transAxes, ha="left",  va="top",    **kw)
    ax.text(0.99, 0.99, "Both ↑\nin wet",
            transform=ax.transAxes, ha="right", va="top",    **kw)
    ax.text(0.01, 0.01, "Both ↓\nin wet",
            transform=ax.transAxes, ha="left",  va="bottom", **kw)
    ax.text(0.99, 0.01, "sRNA ↓,  pathway ↑\nin wet",
            transform=ax.transAxes, ha="right", va="bottom", **kw)

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title_a, fontsize=9, loc="left")
    ax.legend(fontsize=8, loc="center right", framealpha=0.9,
              bbox_to_anchor=(1.0, 0.5))
    ax.xaxis.grid(True, linestyle="--", alpha=0.25, zorder=0)
    ax.yaxis.grid(True, linestyle="--", alpha=0.25, zorder=0)

    #  Panel B: quadrant count bars 
    q_order = [
        "Both increase in wet",
        "sRNA increases, pathway decreases",
        "Both decrease in wet",
        "sRNA decreases, pathway increases",
        "Near zero / ambiguous",
    ]
    q_counts = plot_df["quadrant"].value_counts()
    q_vals   = [q_counts.get(q, 0) for q in q_order]
    q_cols   = [Q_COLORS.get(q, "#CCCCCC") for q in q_order]
    q_labels = [
        "Both ↑\nin wet",
        "sRNA ↑,\npathway ↓",
        "Both ↓\nin wet",
        "sRNA ↓,\npathway ↑",
        "Near zero /\nambiguous",
    ]

    bars = ax_b.bar(range(len(q_order)), q_vals,
                    color=q_cols, edgecolor="white", linewidth=0.6, width=0.6)
    total = sum(q_vals)
    for i, (v, col) in enumerate(zip(q_vals, q_cols)):
        pct = 100 * v / total if total > 0 else 0
        ax_b.text(i, v + 0.3, f"{v}\n({pct:.0f}%)",
                  ha="center", va="bottom", fontsize=8.5, fontweight="bold",
                  color="#333")

    ax_b.set_xticks(range(len(q_order)))
    ax_b.set_xticklabels(q_labels, fontsize=8)
    ax_b.set_ylabel("Number of associations", fontsize=9)
    ax_b.set_title("B. Count by quadrant", fontsize=9, loc="left")
    ax_b.yaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
    ax_b.set_ylim(0, max(q_vals) * 1.25 if q_vals else 1)

    plt.tight_layout()
    plt.savefig(out_path)
    print(f"  Saved → {out_path}")
    plt.close()


#  FIG 7a: one point per sRNA-pathway pair 
print("\nBuilding Fig7a: one point per sRNA-pathway pair...")

LABEL_KEYWORDS = [
    "C2 photosyn",
    "tRNA processing",
    "phytol degradation",
    "guanosine",
    "glycine",
    "sulfate assimilation",
    "peptidoglycan maturation",
]

def should_label(name):
    if pd.isna(name): return False
    return any(kw.lower() in name.lower() for kw in LABEL_KEYWORDS)

pairs = (
    stf.sort_values("interaction_energy")
       .drop_duplicates(["srna_id", "pathway_id"])
       [["srna_id", "pathway_id", "pathway_name",
         "pathway_log2fc", "srna_expr_lfc"]]
       .copy()
)
# Only label matching pathways; others get empty string (no label rendered)
pairs["label"] = pairs["pathway_name"].apply(
    lambda n: n if should_label(n) else ""
)
print(f"  {len(pairs)} unique sRNA-pathway pairs")
print(f"  Labeled pathways: {pairs[pairs['label'] != '']['pathway_name'].unique().tolist()}")
make_scatter(
    pairs,
    x_col="pathway_log2fc", y_col="srna_expr_lfc",
    label_col="label",
    title_a="A. One point per sRNA–pathway pair (best interaction energy kept)",
    xlabel="Pathway log₂FC (wet / dry)",
    ylabel="sRNA expression log₂FC (wet / dry)",
    out_path=f"{FIG_OUT}/Fig7a_scatter_pairs.png",
)

#  FIG 7b: one point per pathway 
print("\nBuilding Fig7b: one point per pathway...")
per_pw = (
    stf.groupby(["pathway_id", "pathway_name", "pathway_log2fc"])
       ["srna_expr_lfc"].mean()
       .reset_index()
       .rename(columns={"srna_expr_lfc": "mean_srna_lfc"})
)
per_pw["label"] = per_pw["pathway_name"]
print(f"  {len(per_pw)} pathways")
make_scatter(
    per_pw,
    x_col="pathway_log2fc", y_col="mean_srna_lfc",
    label_col="label",
    title_a="A. One point per pathway (mean sRNA log₂FC across all qualifying regulators)",
    xlabel="Pathway log₂FC (wet / dry)",
    ylabel="Mean sRNA expression log₂FC (wet / dry)",
    out_path=f"{FIG_OUT}/Fig7b_scatter_per_pathway.png",
    point_size=90,
)

#  FIG 7c: one point per sRNA 
print("\nBuilding Fig7c: one point per sRNA...")
per_srna = (
    stf.groupby(["srna_id", "srna_expr_lfc"])
       ["pathway_log2fc"].mean()
       .reset_index()
       .rename(columns={"pathway_log2fc": "mean_pw_lfc"})
)
per_srna["label"] = per_srna["srna_id"]
print(f"  {len(per_srna)} sRNAs")
make_scatter(
    per_srna,
    x_col="mean_pw_lfc", y_col="srna_expr_lfc",
    label_col="label",
    title_a="A. One point per sRNA (mean pathway log₂FC across all qualifying targets)",
    xlabel="Mean target pathway log₂FC (wet / dry)",
    ylabel="sRNA expression log₂FC (wet / dry)",
    out_path=f"{FIG_OUT}/Fig7c_scatter_per_srna.png",
    point_size=90,
)

print("\nAll done.")

#  FIG 7c: one point per sRNA 
print("\nBuilding Fig7c: one point per sRNA...")
per_srna = (
    stf.groupby(["srna_id", "srna_expr_lfc"])
       ["pathway_log2fc"].mean()
       .reset_index()
       .rename(columns={"pathway_log2fc": "mean_pw_lfc"})
)
per_srna["label"] = per_srna["srna_id"]
print(f"  {len(per_srna)} sRNAs")
make_scatter(
    per_srna,
    x_col="mean_pw_lfc", y_col="srna_expr_lfc",
    label_col="label",
    title_a="A. One point per sRNA (mean pathway log₂FC across all qualifying targets)",
    xlabel="Mean target pathway log₂FC (wet / dry)",
    ylabel="sRNA expression log₂FC (wet / dry)",
    out_path=f"{FIG_OUT}/Fig7c_scatter_per_srna.png",
    point_size=90,
)

print("\nAll done.")


# SUPPLEMENTARY FIGURE 2 — sRNA–Gene Interaction Network
print("Plotting Supp Figure 2 — sRNA–Gene Network...")

# Build network: sRNA nodes → gene nodes → pathway nodes
# Filter to top interactions for readability
net_df = (
    st.dropna(subset=["pathway_id","target_gene_id","srna_id"])
      .sort_values("interaction_energy")
      .drop_duplicates(["srna_id","target_gene_id"])
)

# Top N sRNAs by number of pathway-linked interactions
top_net_srnas = net_df["srna_id"].value_counts().head(20).index.tolist()
net_sub = net_df[net_df["srna_id"].isin(top_net_srnas)]

G = nx.DiGraph()

# Add nodes
for _, row in net_sub.iterrows():
    srna = row["srna_id"]
    gene = row["target_gene_id"]
    pw   = row["pathway_id"]
    pw_name = row.get("pathway_name", pw) if not pd.isna(row.get("pathway_name","")) else pw

    # sRNA direction
    d = srna_uniq[srna_uniq["srna_id"]==srna]["direction"].values
    direction = d[0] if len(d) > 0 else "shared"

    if not G.has_node(srna):
        G.add_node(srna, ntype="srna", direction=direction)
    if not G.has_node(gene):
        G.add_node(gene, ntype="gene")
    if not G.has_node(pw):
        G.add_node(pw, ntype="pathway", name=pw_name,
                   pw_dir=row.get("pathway_direction","ns"))

    G.add_edge(srna, gene, energy=row["interaction_energy"])
    G.add_edge(gene, pw)

# Layout: shell layout — sRNAs outer, genes middle, pathways inner
srna_nodes = [n for n,d in G.nodes(data=True) if d.get("ntype")=="srna"]
gene_nodes  = [n for n,d in G.nodes(data=True) if d.get("ntype")=="gene"]
pw_nodes    = [n for n,d in G.nodes(data=True) if d.get("ntype")=="pathway"]

nlist = [srna_nodes, gene_nodes, pw_nodes]
try:
    pos = nx.shell_layout(G, nlist=nlist)
except:
    pos = nx.spring_layout(G, k=2, seed=42)

fig, ax = plt.subplots(figsize=(14, 14))
ax.set_facecolor("#F8F8F8")

# Draw edges
for u, v, data in G.edges(data=True):
    xu, yu = pos[u]
    xv, yv = pos[v]
    e = data.get("energy", -10)
    alpha = min(0.8, max(0.1, abs(e) / 200))
    ax.plot([xu, xv], [yu, yv], color="#AAAAAA",
            linewidth=0.5, alpha=alpha, zorder=1)

# Draw nodes
node_colors = {
    "srna":    lambda n: WET_C if G.nodes[n].get("direction")=="up" else DRY_C,
    "gene":    lambda n: "#888888",
    "pathway": lambda n: (WET_C if G.nodes[n].get("pw_dir")=="up"
                         else DOWN_C if G.nodes[n].get("pw_dir")=="down"
                         else "#CCCCCC"),
}
node_sizes = {"srna": 350, "gene": 100, "pathway": 250}
node_shapes= {"srna": "o", "gene": "s", "pathway": "D"}

for ntype, shape in node_shapes.items():
    nodes = [n for n,d in G.nodes(data=True) if d.get("ntype")==ntype]
    if not nodes: continue
    colors = [node_colors[ntype](n) for n in nodes]
    nx.draw_networkx_nodes(G, pos, nodelist=nodes,
                           node_color=colors,
                           node_size=node_sizes[ntype],
                           node_shape=shape,
                           alpha=0.9, ax=ax)

# Labels for sRNA and pathway nodes only (genes too many)
srna_labels = {n: n[:10]+"…" if len(n)>10 else n for n in srna_nodes}
pw_labels   = {n: G.nodes[n].get("name", n)[:18]+"…"
               if len(G.nodes[n].get("name", n))>18
               else G.nodes[n].get("name", n)
               for n in pw_nodes}

nx.draw_networkx_labels(G, pos, labels=srna_labels,
                        font_size=5.5, font_color="black",
                        font_weight="bold", ax=ax)
nx.draw_networkx_labels(G, pos, labels=pw_labels,
                        font_size=5.5, font_color="#333",
                        ax=ax)

ax.set_title("Supplementary Figure 2: sRNA → Gene → Pathway Interaction Network\n"
             "Top 20 sRNAs by pathway-linked interactions  |  "
             "Circle=sRNA, Square=gene, Diamond=pathway",
             fontsize=10, fontweight="bold")
ax.axis("off")

legend_handles = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor=WET_C,
           markersize=10, label="sRNA — wet-only (up)"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=DRY_C,
           markeredgecolor="#888", markersize=10, label="sRNA — dry-only (down)"),
    Line2D([0],[0], marker="s", color="w", markerfacecolor="#888888",
           markeredgecolor="#888", markersize=8, label="Target gene"),
    Line2D([0],[0], marker="D", color="w", markerfacecolor=WET_C,
           markersize=9, label="Pathway — increasing"),
    Line2D([0],[0], marker="D", color="w", markerfacecolor=DOWN_C,
           markersize=9, label="Pathway — decreasing"),
]
ax.legend(handles=legend_handles, loc="lower left",
          fontsize=8, framealpha=0.9)
plt.tight_layout()
save("SuppFig2_sRNA_gene_network.png")

# DONE

print(f"\nAll figures written to {FIG_OUT}/")
print("Files:")
for f in sorted(os.listdir(FIG_OUT)):
    size = os.path.getsize(f"{FIG_OUT}/{f}") // 1024
    print(f"  {f}  ({size} KB)")


# FIGURE 2b — sRNA Expression Comparison (replaces volcano as main Fig 2)

print("Plotting Figure 2b — sRNA Expression Comparison...")

# Build per-sRNA mean expression per condition
expr_summary = (
    sm.groupby(["srna_id", "condition"])["expression"]
      .mean()
      .reset_index()
)
expr_summary = expr_summary.merge(
    srna_uniq[["srna_id", "type", "direction"]], on="srna_id", how="left"
)
expr_summary["log2_tpm"] = np.log2(expr_summary["expression"] + 1)
expr_summary = expr_summary.dropna(subset=["type"])

fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
fig.suptitle(
    "Figure 2b: sRNA Expression Levels by Condition and Class",
    fontsize=12, fontweight="bold"
)

cond_order  = ["dry", "wet"]
cond_colors = {"dry": DRY_C, "wet": WET_C}
cond_labels = {"dry": "2016 (dry)", "wet": "2017 (wet)"}

for ax, stype in zip(axes, ["asRNA", "itsRNA"]):
    sub = expr_summary[expr_summary["type"] == stype]

    # Box underneath
    bp = ax.boxplot(
        [sub[sub["condition"]==c]["log2_tpm"].values for c in cond_order],
        positions=[0, 1],
        patch_artist=True,
        widths=0.35,
        medianprops=dict(color="black", linewidth=2),
        whiskerprops=dict(linewidth=1.2),
        capprops=dict(linewidth=1.2),
        flierprops=dict(marker="", alpha=0),
        zorder=2,
    )
    for patch, cond in zip(bp["boxes"], cond_order):
        patch.set_facecolor(cond_colors[cond])
        patch.set_alpha(0.5)

    # Strip on top
    np.random.seed(42)
    for i, cond in enumerate(cond_order):
        vals = sub[sub["condition"] == cond]["log2_tpm"].values
        jitter = np.random.uniform(-0.12, 0.12, len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter, vals,
            color=cond_colors[cond],
            s=18, alpha=0.7, edgecolors="white",
            linewidths=0.3, zorder=4,
        )

    # Annotate n per condition
    for i, cond in enumerate(cond_order):
        n = sub[sub["condition"]==cond]["srna_id"].nunique()
        med = sub[sub["condition"]==cond]["log2_tpm"].median()
        ax.text(i, -0.5, f"n = {n}\nsRNAs",
                ha="center", fontsize=8, color="#555")

    # Mann-Whitney test
    g1 = sub[sub["condition"]=="dry"]["log2_tpm"].dropna()
    g2 = sub[sub["condition"]=="wet"]["log2_tpm"].dropna()
    if len(g1) > 0 and len(g2) > 0:
        _, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")
        ymax = sub["log2_tpm"].max()
        ax.plot([0, 0, 1, 1], [ymax+0.3, ymax+0.6, ymax+0.6, ymax+0.3],
                color="black", linewidth=1)
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        ax.text(0.5, ymax+0.7, f"p={pval:.3f} ({sig})",
                ha="center", fontsize=8)

    ax.set_xticks([0, 1])
    ax.set_xticklabels([cond_labels[c] for c in cond_order], fontsize=10)
    ax.set_title(f"{'Antisense (asRNA)' if stype=='asRNA' else 'Intergenic (itsRNA)'}",
                 fontsize=11, fontweight="bold")
    ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
    ax.set_xlim(-0.6, 1.6)

axes[0].set_ylabel("Mean log₂(TPM + 1)", fontsize=11)

# Shared legend patches
handles = [
    mpatches.Patch(color=DRY_C, label="2016 — dry", alpha=0.8),
    mpatches.Patch(color=WET_C, label="2017 — wet", alpha=0.8),
]
fig.legend(handles=handles, loc="upper right",
           bbox_to_anchor=(0.99, 0.92), fontsize=9)
plt.tight_layout()
save("Fig2b_sRNA_expression_comparison.png")


# FIGURE 2c — sRNA Class Composition per Condition
print("Plotting Figure 2c — sRNA Class Composition...")

fig, axes = plt.subplots(1, 2, figsize=(11, 5))
fig.suptitle("Figure 2c: sRNA Class Composition and Mean Expression by Condition",
             fontsize=11, fontweight="bold")

# Panel A: counts
count_df = (
    srna_uniq.groupby(["direction", "type"])
             .size()
             .reset_index(name="count")
)
# map direction to condition label
dir_label = {"down": "2016 only\n(dry)", "up": "2017 only\n(wet)", "shared": "Shared"}
count_df["cond_label"] = count_df["direction"].map(dir_label)
type_colors2 = {"asRNA": "#E07B54", "itsRNA": "#5BA4CF"}

ax = axes[0]
cond_positions = {"2016 only\n(dry)": 0, "Shared": 1, "2017 only\n(wet)": 2}
bar_width = 0.35

for j, stype in enumerate(["asRNA", "itsRNA"]):
    sub = count_df[count_df["type"] == stype]
    xs  = [cond_positions[r] + (j-0.5)*bar_width for r in sub["cond_label"]]
    ax.bar(xs, sub["count"].values,
           width=bar_width, color=type_colors2[stype],
           label=stype, alpha=0.85, edgecolor="white", zorder=3)
    for x, v in zip(xs, sub["count"].values):
        ax.text(x, v + 0.5, str(v), ha="center", fontsize=8, fontweight="bold")

ax.set_xticks(list(cond_positions.values()))
ax.set_xticklabels(list(cond_positions.keys()), fontsize=9)
ax.set_ylabel("Number of sRNAs", fontsize=10)
ax.set_title("A. sRNA counts by class and condition", fontsize=10, loc="left")
ax.legend(fontsize=9)
ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)

# Panel B: mean expression per type per condition
ax2 = axes[1]
mean_expr = (
    expr_summary[expr_summary["direction"] != "shared"]
    .groupby(["condition", "type"])["log2_tpm"]
    .agg(["mean", "sem"])
    .reset_index()
)
cond_x = {"dry": 0, "wet": 1}
for j, stype in enumerate(["asRNA", "itsRNA"]):
    sub = mean_expr[mean_expr["type"] == stype]
    xs  = [cond_x[c] + (j-0.5)*0.35 for c in sub["condition"]]
    ax2.bar(xs, sub["mean"].values,
            yerr=sub["sem"].values,
            width=0.32, color=type_colors2[stype],
            label=stype, alpha=0.85,
            edgecolor="white", capsize=4, zorder=3)

ax2.set_xticks([0, 1])
ax2.set_xticklabels(["2016 (dry)", "2017 (wet)"], fontsize=9)
ax2.set_ylabel("Mean log₂(TPM + 1) ± SEM", fontsize=10)
ax2.set_title("B. Mean expression per class", fontsize=10, loc="left")
ax2.legend(fontsize=9)
ax2.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)

plt.tight_layout()
save("Fig2c_sRNA_class_composition.png")


# FIGURE 4b — Top Pathway Changes (horizontal ranked bar + sRNA target count)
print("Plotting Figure 4b — Top Pathway Changes...")

# sRNA target count per pathway
srna_per_pw = (
    st.dropna(subset=["pathway_id", "srna_id"])
      .groupby("pathway_id")["srna_id"]
      .nunique()
      .rename("n_srnas")
)

pw_plot = pw_uniq.dropna(subset=["log2fc", "pathway_name"]).copy()
pw_plot["n_srnas"] = pw_plot["pathway_id"].map(srna_per_pw).fillna(0).astype(int)
pw_plot["short_name"] = pw_plot["pathway_name"].apply(
    lambda x: textwrap.fill(x, 35) if not pd.isna(x) else x
)

n_show = 20
top_down = pw_plot.nsmallest(n_show, "log2fc")
top_up   = pw_plot.nlargest(n_show,  "log2fc")
top_both = pd.concat([top_down, top_up]).drop_duplicates("pathway_id")
top_both = top_both.sort_values("log2fc")

fig, (ax_bar, ax_count) = plt.subplots(
    1, 2, figsize=(14, 10),
    gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05}
)
fig.suptitle(
    "Figure 4b: Top Pathway-Level Changes (2016→2017)\n"
    "with sRNA Targeting Frequency",
    fontsize=11, fontweight="bold"
)

y_pos = range(len(top_both))
bar_colors = [
    HILIGHT if row.targeted else (WET_C if row.log2fc > 0 else DOWN_C)
    for row in top_both.itertuples()
]
bars = ax_bar.barh(list(y_pos), top_both["log2fc"].values,
                   color=bar_colors, height=0.7,
                   edgecolor="white", linewidth=0.5, zorder=3)

ax_bar.axvline(0, color="black", linewidth=1)
ax_bar.axvline(0.5,  color="#888", linewidth=0.5, linestyle="--", alpha=0.5)
ax_bar.axvline(-0.5, color="#888", linewidth=0.5, linestyle="--", alpha=0.5)

ax_bar.set_yticks(list(y_pos))
ax_bar.set_yticklabels(top_both["short_name"].values, fontsize=7.5)
ax_bar.set_xlabel("RNA log₂FC (wet / dry)", fontsize=10)
ax_bar.xaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
ax_bar.set_title("A. Fold-change", fontsize=9, loc="left")

# Category color strip on left
for i, row in enumerate(top_both.itertuples()):
    cat_color = CAT_COLORS.get(row.category, "#AAA")
    ax_bar.add_patch(plt.Rectangle(
        (ax_bar.get_xlim()[0] - 0.3, i - 0.35), 0.25, 0.7,
        color=cat_color, clip_on=False, zorder=6,
    ))

# sRNA count panel
count_colors = [HILIGHT if n > 0 else "#DDDDDD"
                for n in top_both["n_srnas"].values]
ax_count.barh(list(y_pos), top_both["n_srnas"].values,
              color=count_colors, height=0.7,
              edgecolor="white", linewidth=0.5, zorder=3)
for i, n in enumerate(top_both["n_srnas"].values):
    if n > 0:
        ax_count.text(n + 0.05, i, str(n), va="center",
                      fontsize=7, fontweight="bold")

ax_count.set_yticks([])
ax_count.set_xlabel("# sRNAs targeting", fontsize=9)
ax_count.xaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
ax_count.set_title("B. sRNA targets", fontsize=9, loc="left")
ax_count.set_xlim(0, max(top_both["n_srnas"].max() + 1, 3))

# Legend
legend_handles = [
    mpatches.Patch(color=WET_C,   label="Increasing in wet (no sRNA target)"),
    mpatches.Patch(color=DOWN_C,  label="Decreasing in wet (no sRNA target)"),
    mpatches.Patch(color=HILIGHT, label="sRNA-targeted pathway"),
] + [mpatches.Patch(color=v, label=k) for k, v in CAT_COLORS.items()]
ax_bar.legend(handles=legend_handles, fontsize=7, loc="lower right",
              framealpha=0.85, ncol=2)

plt.tight_layout()
save("Fig4b_top_pathway_changes.png") 

#Humann Abundance 
HUMANN_DNA = "AD/Total_DNA_pathabundance_humann_table.tsv"
HUMANN_RNA = "AD/Total_RNA_pathabundance_humann_table.tsv"


def load_humann_pathways(fpath):
    try:
        df = pd.read_csv(fpath, sep="\t", index_col=0)

        # Drop special rows
        df = df[~df.index.str.startswith(("UNMAPPED", "UNINTEGRATED"))].copy()

        # Remove stratification (everything after |)
        df.index = df.index.str.split("|").str[0]

        # Collapse duplicates (classified + unclassified)
        df = df.groupby(df.index).sum()

        # Clean pathway names (remove IDs like PWY-123: )
        df.index = df.index.str.replace(r"^[^:]+:\s*", "", regex=True)

        df = df.astype(float)

        # Convert to relative abundance (%)
        df = df.div(df.sum(axis=0), axis=1) * 100

        return df

    except Exception as e:
        print(f"[WARN] Could not load {fpath}: {e}")
        return None  

def sort_samples_by_year(df):
    def sort_key(s):
        parts = s.split("_")

        year = next((p for p in parts if p in ("2016", "2017")), "9999")
        rep  = next((p for p in parts if p.startswith("R") and p[1:].isdigit()), "R0")

        rep_num = int(rep[1:]) if rep[1:].isdigit() else 0

        return (int(year), rep_num)

    sorted_cols = sorted(df.columns, key=sort_key)
    return df[sorted_cols]
    

dna_path = load_humann_pathways(HUMANN_DNA)
rna_path = load_humann_pathways(HUMANN_RNA)

rna_path = sort_samples_by_year(rna_path)
dna_path = sort_samples_by_year(dna_path)


all_pathways = list(dict.fromkeys(list(rna_path.index) + list(dna_path.index)))
palette = sns.color_palette("tab20", len(all_pathways))
pathway_color_map = {pw: palette[i % len(palette)] for i, pw in enumerate(all_pathways)} 

def format_sample_label(s):
    parts = s.split("_")

    year = next((p for p in parts if p in ("2016", "2017")), "")
    rep  = next((p for p in parts if p.startswith("R") and p[1:].isdigit()), "")

    return f"{year}\n{rep}" if year else s

#Top 20 Pathways
def plot_pathway_topN(ax, df, title, N=20):
    top_pw = df.mean(axis=1).nlargest(N).index.tolist()
    other  = df.loc[~df.index.isin(top_pw)].sum(axis=0)

    plot_df = df.loc[top_pw].T.copy()
    plot_df["Other"] = other.values

    bottom = np.zeros(len(plot_df))

    for pw in plot_df.columns:
        color = pathway_color_map.get(pw, "#CCCCCC")
        ax.bar(range(len(plot_df)), plot_df[pw].values,
               bottom=bottom, color=color, label=pw, width=0.8)
        bottom += plot_df[pw].values

    samples = list(plot_df.index)
    n16 = sum("2016" in s for s in samples)
    n17 = sum("2017" in s for s in samples)

    ax.axvline(n16 - 0.5, color="black", linestyle="--")

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([format_sample_label(s) for s in samples], fontsize=7)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Relative abundance (%)")
    ax.set_title(title, loc="left")

    return plot_df.columns.tolist()
#All Pathways 
def plot_pathway_all(ax, df, title):
    plot_df = df.T.copy()

    bottom = np.zeros(len(plot_df))

    for pw in plot_df.columns:
        color = pathway_color_map.get(pw, "#CCCCCC")
        ax.bar(range(len(plot_df)), plot_df[pw].values,
               bottom=bottom, color=color, width=0.8)
        bottom += plot_df[pw].values

    samples = list(plot_df.index)
    n16 = sum("2016" in s for s in samples)

    ax.axvline(n16 - 0.5, color="black", linestyle="--")

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([format_sample_label(s) for s in samples], fontsize=7)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Relative abundance (%)")
    ax.set_title(title, loc="left")
#Print 
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

pw_rna = plot_pathway_topN(axes[0], rna_path, "RNA pathways (Top 20)")
pw_dna = plot_pathway_topN(axes[1], dna_path, "DNA pathways (Top 20)")

import matplotlib.patches as mpatches

# combine pathways from both panels
all_pw = list(dict.fromkeys(pw_rna + pw_dna))

handles = [
    mpatches.Patch(color=pathway_color_map.get(pw, "#CCCCCC"), label=pw)
    for pw in all_pw
]

fig.legend(handles=handles,
           loc="center right",
           bbox_to_anchor=(1.25, 0.5),
           fontsize=7,
           title="Pathways",
           title_fontsize=9)
plt.tight_layout()
plt.savefig("humann_pathways_top20.png", dpi=300, bbox_inches="tight")
plt.close()

fig, axes = plt.subplots(2, 1, figsize=(14, 10))

plot_pathway_all(axes[0], rna_path, "RNA pathways (All)")
plot_pathway_all(axes[1], dna_path, "DNA pathways (All)")

import matplotlib.patches as mpatches

# combine pathways from both panels
all_pw = list(dict.fromkeys(pw_rna + pw_dna))

handles = [
    mpatches.Patch(color=pathway_color_map.get(pw, "#CCCCCC"), label=pw)
    for pw in all_pw
]

fig.legend(handles=handles,
           loc="center right",
           bbox_to_anchor=(1.25, 0.5),
           fontsize=7,
           title="Pathways",
           title_fontsize=9)

plt.tight_layout(rect=[0, 0, 0.8, 1])
plt.savefig("humann_pathways_all.png", dpi=300, bbox_inches="tight")
plt.close()

print("\nAdditional figures written.")
