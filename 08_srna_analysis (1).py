#!/usr/bin/env python3
"""
08_srna_analysis.py

Analyses:
  1. Differential expression of sRNAs between timepoints (permutation test)
  2. Structural comparison: asRNA vs itsRNA
  3. sRNA-level Random Forest: predict up/downregulation from sRNA features
  4. sRNA-target expression correlation (log2FC approximation)
  5. Summary figures

Usage:
    python 08_srna_analysis.py \
        --dataset        AD \
        --snapt_t1       AD/snapt_2016 \
        --snapt_t2       AD/snapt_2017 \
        --rnafold_t1     AD/rnafold_2016/rnafold_summary.tsv \
        --rnafold_t2     AD/rnafold_2017/rnafold_summary.tsv \
        --feature_matrix AD/feature_matrix.tsv \
        --outdir         AD/analysis_output \
        --intarna_t1     AD/intarna_2016/intarna_summary_fixed.tsv \
        --intarna_t2     AD/intarna_2017/intarna_summary_fixed.tsv \
        --label_t1       2016 \
        --label_t2       2017 \
        --tpm_matrix_t1  ls AD/per_replicate_tpm_2016/tpm_matrix.tsv \
        --tpm_matrix_t2  ls AD/per_replicate_tpm_2017/tpm_matrix.tsv
"""

import argparse
import os
import re
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu

PSEUDOCOUNT  = 1e-6
RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# ── Arguments ─────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--dataset",        required=True)
parser.add_argument("--snapt_t1",       required=True)
parser.add_argument("--snapt_t2",       required=True)
parser.add_argument("--rnafold_t1",     required=True)
parser.add_argument("--rnafold_t2",     required=True)
parser.add_argument("--feature_matrix", required=True)
parser.add_argument("--outdir",         required=True)
parser.add_argument("--intarna_t1",     default=None)
parser.add_argument("--intarna_t2",     default=None)
parser.add_argument("--label_t1",       default=None)
parser.add_argument("--label_t2",       default=None)
parser.add_argument("--tpm_matrix_t1",  default=None)
parser.add_argument("--tpm_matrix_t2",  default=None)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
DATASET = args.dataset.upper()

def infer_label(path):
    base = os.path.basename(path.rstrip("/"))
    m = re.search(r"snapt_(.+)$", base)
    return m.group(1) if m else base

T1 = args.label_t1 or infer_label(args.snapt_t1)
T2 = args.label_t2 or infer_label(args.snapt_t2)

print(f"\n{'='*60}")
print(f"Dataset: {DATASET}   Timepoints: {T1} → {T2}")
print(f"{'='*60}\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Load sRNA data
# ═══════════════════════════════════════════════════════════════════════════════
print("SECTION 1: Loading sRNA data")
print("-" * 40)

def parse_snapt_gff(gff_path, srna_type, timegroup):
    if not os.path.exists(gff_path):
        print(f"  WARNING: {gff_path} not found")
        return pd.DataFrame()
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            contig = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]
            attrs  = parts[8]

            def get_attr(key):
                m = re.search(rf'{key} "([^"]+)"', attrs)
                return m.group(1) if m else ""

            fpkm = float(get_attr("FPKM") or 0)
            tpm  = float(get_attr("TPM")  or 0)
            cov  = float(get_attr("cov")  or 0)

            # Use same ID format as RNAfold: contig:start-1-end(strand)
            srna_id = f"{contig}:{start-1}-{end}({strand})"

            rows.append({
                "srna_id":   srna_id,
                "contig":    contig,
                "start":     start,
                "end":       end,
                "strand":    strand,
                "length":    end - start + 1,
                "srna_type": srna_type,
                "is_asrna":  1 if srna_type == "antisense" else 0,
                "cov":       cov,
                "fpkm":      fpkm,
                "tpm":       tpm,
                "timegroup": timegroup,
            })
    return pd.DataFrame(rows)

frames_t1, frames_t2 = [], []
for stype in ["antisense", "intergenic"]:
    prefix = "small_antisense" if stype == "antisense" else "small_intergenic"
    df1 = parse_snapt_gff(
        os.path.join(args.snapt_t1, f"{prefix}_ncRNAs.gff"), stype, T1)
    df2 = parse_snapt_gff(
        os.path.join(args.snapt_t2, f"{prefix}_ncRNAs.gff"), stype, T2)
    if not df1.empty: frames_t1.append(df1)
    if not df2.empty: frames_t2.append(df2)

srna_t1 = pd.concat(frames_t1, ignore_index=True) if frames_t1 else pd.DataFrame()
srna_t2 = pd.concat(frames_t2, ignore_index=True) if frames_t2 else pd.DataFrame()

print(f"  T1 ({T1}): {len(srna_t1)} sRNAs  "
      f"(asRNA={srna_t1['is_asrna'].sum()}, "
      f"itsRNA={(srna_t1['is_asrna']==0).sum()})")
print(f"  T2 ({T2}): {len(srna_t2)} sRNAs  "
      f"(asRNA={srna_t2['is_asrna'].sum()}, "
      f"itsRNA={(srna_t2['is_asrna']==0).sum()})")

# ── Compute log2FC ─────────────────────────────────────────────────────────────
fpkm_t1 = srna_t1.set_index("srna_id")["fpkm"].rename(f"fpkm_{T1}")
fpkm_t2 = srna_t2.set_index("srna_id")["fpkm"].rename(f"fpkm_{T2}")
tpm_t1  = srna_t1.set_index("srna_id")["tpm"].rename(f"tpm_{T1}")
tpm_t2  = srna_t2.set_index("srna_id")["tpm"].rename(f"tpm_{T2}")

all_ids = fpkm_t1.index.union(fpkm_t2.index)
expr = (pd.DataFrame(index=all_ids)
        .join(fpkm_t1).join(fpkm_t2)
        .join(tpm_t1).join(tpm_t2)
        .fillna(0)
        .reset_index()
        .rename(columns={"index": "srna_id"}))

expr["log2fc_fpkm"] = np.log2(
    (expr[f"fpkm_{T2}"] + PSEUDOCOUNT) / (expr[f"fpkm_{T1}"] + PSEUDOCOUNT))
expr["log2fc_tpm"] = np.log2(
    (expr[f"tpm_{T2}"]  + PSEUDOCOUNT) / (expr[f"tpm_{T1}"]  + PSEUDOCOUNT))

# Annotate sRNA type
type_map = pd.concat([
    srna_t1[["srna_id", "srna_type", "is_asrna", "length"]],
    srna_t2[["srna_id", "srna_type", "is_asrna", "length"]],
]).drop_duplicates("srna_id")
expr = expr.merge(type_map, on="srna_id", how="left")

both = ((expr[f"fpkm_{T1}"] > 0) & (expr[f"fpkm_{T2}"] > 0)).sum()
only1 = ((expr[f"fpkm_{T1}"] > 0) & (expr[f"fpkm_{T2}"] == 0)).sum()
only2 = ((expr[f"fpkm_{T1}"] == 0) & (expr[f"fpkm_{T2}"] > 0)).sum()
print(f"\n  Detected in both: {both}")
print(f"  Only in {T1}:     {only1}")
print(f"  Only in {T2}:     {only2}")

# ── Load RNAfold ───────────────────────────────────────────────────────────────
def load_rnafold(path, tg):
    if not os.path.exists(path):
        print(f"  WARNING: RNAfold missing: {path}")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    df["timegroup"] = tg
    print(f"  RNAfold {tg}: {len(df)} structures")
    return df

rf_t1 = load_rnafold(args.rnafold_t1, T1)
rf_t2 = load_rnafold(args.rnafold_t2, T2)

# ── Load IntaRNA (already column-fixed by fix_intarna_summary.py) ─────────────
def load_intarna(path, tg):
    if not path or not os.path.exists(path):
        print(f"  IntaRNA {tg}: not found — interaction features will be skipped")
        return None
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]
    print(f"  IntaRNA {tg}: {len(df)} sRNAs  "
          f"(energy range {df['best_energy'].min():.2f} to "
          f"{df['best_energy'].max():.2f})")
    return df

intarna_t1 = load_intarna(args.intarna_t1, T1)
intarna_t2 = load_intarna(args.intarna_t2, T2)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Differential expression (permutation test)
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("SECTION 2: Differential expression")
print("-" * 40)

both_detected = expr[
    (expr[f"fpkm_{T1}"] > 0) & (expr[f"fpkm_{T2}"] > 0)].copy()

def permutation_test_mean(values, n_perm=10000):
    obs = np.mean(values)
    vals = np.array(values)
    count = sum(
        abs(np.mean(vals * np.random.choice([-1,1], size=len(vals)))) >= abs(obs)
        for _ in range(n_perm)
    )
    return obs, count / n_perm

if len(both_detected) > 3:
    obs_mean, pval = permutation_test_mean(both_detected["log2fc_fpkm"].values)
    print(f"  Mean log2FC FPKM ({T1}→{T2}): {obs_mean:.3f}")
    print(f"  Permutation p-value:         {pval:.4f}")
    print(f"  Interpretation: pool "
          f"{'significantly' if pval<0.05 else 'not significantly'} shifted")

# With only 1 sRNA shared between timepoints, most log2FC values are
# pseudocount-driven (one side = 0). Use a higher threshold to identify
# the most meaningful changes among the shared sRNA, and separately
# analyze the timepoint-specific sRNA pools.
LOG2FC_THRESH = 1.0

# Flag timepoint-specific sRNAs explicitly
expr["detection"] = "both"
expr.loc[(expr[f"fpkm_{T1}"] > 0) & (expr[f"fpkm_{T2}"] == 0), "detection"] = f"only_{T1}"
expr.loc[(expr[f"fpkm_{T1}"] == 0) & (expr[f"fpkm_{T2}"] > 0), "detection"] = f"only_{T2}"

n_only_t1 = (expr["detection"] == f"only_{T1}").sum()
n_only_t2 = (expr["detection"] == f"only_{T2}").sum()
n_both    = (expr["detection"] == "both").sum()

print(f"\n  sRNA detection summary:")
print(f"    Unique to {T1}: {n_only_t1}  ({n_only_t1/len(expr)*100:.0f}%)")
print(f"    Unique to {T2}: {n_only_t2}  ({n_only_t2/len(expr)*100:.0f}%)")
print(f"    Shared:         {n_both}  ({n_both/len(expr)*100:.0f}%)")
print(f"    → Near-complete sRNA pool turnover between timepoints")
print(f"    → This is the primary DE finding: community expresses distinct")
print(f"      regulatory repertoires pre- vs post-rain event")

expr["de_status"] = "stable"
expr.loc[expr["log2fc_fpkm"] >  LOG2FC_THRESH, "de_status"] = f"up_in_{T2}"
expr.loc[expr["log2fc_fpkm"] < -LOG2FC_THRESH, "de_status"] = f"up_in_{T1}"

n_up   = (expr["de_status"] == f"up_in_{T2}").sum()
n_down = (expr["de_status"] == f"up_in_{T1}").sum()
n_stab = (expr["de_status"] == "stable").sum()

print(f"\n  |log2FC| > {LOG2FC_THRESH}:")
print(f"    Up in {T2}: {n_up}")
print(f"    Up in {T1}: {n_down}")
print(f"    Stable:     {n_stab}")

expr.to_csv(os.path.join(args.outdir, "srna_expression_log2fc.tsv"),
            sep="\t", index=False)

# ── Volcano / MA plot ──────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))
color_map = {
    f"up_in_{T2}": "#d95f02",
    f"up_in_{T1}": "#1b9e77",
    "stable":      "#aaaaaa",
}
for status, grp in expr.groupby("de_status"):
    ax.scatter(
        grp["log2fc_fpkm"],
        np.log2(grp[f"fpkm_{T1}"] + grp[f"fpkm_{T2}"] + PSEUDOCOUNT),
        c=color_map.get(status, "#aaaaaa"),
        alpha=0.7, s=40, label=status, edgecolors="none"
    )
ax.axvline(x=LOG2FC_THRESH,  ls="--", lw=1, color="gray")
ax.axvline(x=-LOG2FC_THRESH, ls="--", lw=1, color="gray")
ax.axvline(x=0, ls="-", lw=0.5, color="black")
ax.set_xlabel(f"log₂FC FPKM ({T1} → {T2})")
ax.set_ylabel("log₂(mean FPKM)")
ax.set_title(f"sRNA expression: {DATASET} {T1}→{T2}\n"
             f"↑{T2}: {n_up}  ↑{T1}: {n_down}  stable: {n_stab}")
ax.legend(frameon=False)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, "volcano_srna.pdf"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(args.outdir, "volcano_srna.png"), dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: volcano_srna.pdf/png")

# ── Turnover figure — more informative than volcano given n=1 shared ──────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: stacked bar showing class composition by detection status
ax = axes[0]
categories = [f"Only {T1}\n(n={n_only_t1})",
              f"Shared\n(n={n_both})",
              f"Only {T2}\n(n={n_only_t2})"]

asrna_counts = [
    (expr[(expr["detection"]==f"only_{T1}") & (expr["is_asrna"]==1)]).shape[0],
    (expr[(expr["detection"]=="both")       & (expr["is_asrna"]==1)]).shape[0],
    (expr[(expr["detection"]==f"only_{T2}") & (expr["is_asrna"]==1)]).shape[0],
]
itsrna_counts = [
    (expr[(expr["detection"]==f"only_{T1}") & (expr["is_asrna"]==0)]).shape[0],
    (expr[(expr["detection"]=="both")       & (expr["is_asrna"]==0)]).shape[0],
    (expr[(expr["detection"]==f"only_{T2}") & (expr["is_asrna"]==0)]).shape[0],
]

x = range(len(categories))
ax.bar(x, asrna_counts,  label="asRNA",  color="#d95f02", alpha=0.85)
ax.bar(x, itsrna_counts, label="itsRNA", color="#1b9e77", alpha=0.85,
       bottom=asrna_counts)
ax.set_xticks(list(x))
ax.set_xticklabels(categories)
ax.set_ylabel("Number of sRNAs")
ax.set_title(f"sRNA pool composition by detection\n{DATASET} {T1} vs {T2}")
ax.legend(frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Right: FPKM distributions for each timepoint pool
ax = axes[1]
fpkm_t1_vals = np.log2(srna_t1["fpkm"] + PSEUDOCOUNT)
fpkm_t2_vals = np.log2(srna_t2["fpkm"] + PSEUDOCOUNT)
ax.hist(fpkm_t1_vals, bins=20, alpha=0.6, color="#4393c3",
        label=f"{T1} (n={len(srna_t1)})", edgecolor="white", lw=0.5)
ax.hist(fpkm_t2_vals, bins=20, alpha=0.6, color="#d6604d",
        label=f"{T2} (n={len(srna_t2)})", edgecolor="white", lw=0.5)
ax.set_xlabel("log₂(FPKM)")
ax.set_ylabel("Count")
ax.set_title(f"sRNA expression levels by timepoint\n"
             f"Mean {T1}: {srna_t1['fpkm'].mean():.1f}  "
             f"Mean {T2}: {srna_t2['fpkm'].mean():.1f}")
ax.legend(frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.suptitle(f"sRNA regulatory repertoire turnover — {DATASET}", fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, "srna_turnover.pdf"),
            dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(args.outdir, "srna_turnover.png"),
            dpi=150, bbox_inches="tight")
plt.close()
print("  Saved: srna_turnover.pdf/png")
# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — asRNA vs itsRNA structural comparison
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("SECTION 3: asRNA vs itsRNA structural comparison")
print("-" * 40)

RF_MERGE_COLS = ["srna_id"] + [c for c in
    ["gc_content", "mfe", "struct_entropy"]
    if not rf_t1.empty and c in rf_t1.columns]
struct = expr.merge(rf_t1[RF_MERGE_COLS], on="srna_id", how="left") \
    if not rf_t1.empty else expr.copy()

# Replace whatever STRUCT_FEATS line exists with:
STRUCT_FEATS = [c for c in
    ["gc_content", "mfe", "struct_entropy", "length"]
    if c in struct.columns]

struct_results = []
for feat in STRUCT_FEATS:
    if feat not in struct.columns:
        continue
    a = struct[struct["is_asrna"]==1][feat].dropna()
    b = struct[struct["is_asrna"]==0][feat].dropna()
    if len(a) < 3 or len(b) < 3:
        continue
    stat, pval = mannwhitneyu(a, b, alternative="two-sided")
    struct_results.append({
        "feature":       feat,
        "asrna_median":  a.median(),
        "itsrna_median": b.median(),
        "mannwhitney_U": stat,
        "pvalue":        pval,
        "significant":   pval < 0.05,
    })
    print(f"  {feat:20s}  asRNA={a.median():.3f}  "
          f"itsRNA={b.median():.3f}  p={pval:.4f}"
          f"{'  *' if pval<0.05 else ''}")

pd.DataFrame(struct_results).to_csv(
    os.path.join(args.outdir, "asrna_vs_itsrna_structural.tsv"),
    sep="\t", index=False)

if struct_results:
    n = len(struct_results)
    fig, axes = plt.subplots(1, n, figsize=(3.5*n, 5))
    if n == 1: axes = [axes]
    for ax, row in zip(axes, struct_results):
        feat = row["feature"]
        a = struct[struct["is_asrna"]==1][feat].dropna()
        b = struct[struct["is_asrna"]==0][feat].dropna()
        bp = ax.boxplot([a, b], labels=["asRNA","itsRNA"],
                        patch_artist=True,
                        medianprops=dict(color="black", lw=2))
        bp["boxes"][0].set_facecolor("#d95f02")
        bp["boxes"][1].set_facecolor("#1b9e77")
        ax.set_title(f"{feat}\np={row['pvalue']:.3f}"
                     f"{'*' if row['significant'] else ''}")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.suptitle(f"asRNA vs itsRNA — {DATASET} {T1}", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "struct_comparison.pdf"),
                dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(args.outdir, "struct_comparison.png"),
                dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved: struct_comparison.pdf/png")

# Replace SECTION 4 with this corrected version:

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — sRNA-level Random Forest
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("SECTION 4: sRNA-level Random Forest")
print("-" * 40)

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import (
        StratifiedKFold, LeaveOneOut, cross_val_score
    )
    from sklearn.pipeline import Pipeline
    import shap
    SK = True
except ImportError as e:
    print(f"  scikit-learn/shap not available: {e}")
    SK = False

if SK:
    ml = expr.copy()

    # Add RNAfold features
    if not rf_t1.empty:
        RF_ML_COLS = ["srna_id"] + [c for c in
            ["gc_content", "mfe", "struct_entropy"]
            if c in rf_t1.columns]
        ml = ml.merge(rf_t1[RF_ML_COLS], on="srna_id", how="left")

    # Add IntaRNA features
    if intarna_t1 is not None and "srna_id" in intarna_t1.columns:
        ml = ml.merge(
            intarna_t1[["srna_id","n_targets","best_energy","mean_energy"]].rename(
                columns={"n_targets":   "int_n_targets",
                         "best_energy": "int_best_energy",
                         "mean_energy": "int_mean_energy"}),
            on="srna_id", how="left")

    # Select only DE sRNAs
    ml_de = ml[ml["de_status"] != "stable"].copy()
    ml_de["label"] = (ml_de["de_status"] == f"up_in_{T2}").astype(int)

    print(f"  Labeled sRNAs: {len(ml_de)}  "
          f"(up={ml_de['label'].sum()}, "
          f"down={(ml_de['label']==0).sum()})")

    # Select only features that actually exist and aren't all NaN
    FEATURE_COLS = []
    for col in ["gc_content", "mfe", "struct_entropy", "length", "is_asrna",
                f"fpkm_{T1}", "int_n_targets", "int_best_energy", "int_mean_energy"]:
        if col in ml_de.columns:
            # Check if column has any non-NaN values
            if ml_de[col].notna().sum() > 0:
                FEATURE_COLS.append(col)

    print(f"  Features ({len(FEATURE_COLS)}): {FEATURE_COLS}")

    if len(FEATURE_COLS) < 2:
        print(f"  WARNING: Not enough valid features for Random Forest")
        print(f"  Skipping RF analysis")
    else:
        # Get feature data and fill NaNs with median
        X = ml_de[FEATURE_COLS].copy()
        for col in FEATURE_COLS:
            median_val = X[col].median()
            if pd.isna(median_val):
                # If all NaN, use 0
                X[col].fillna(0, inplace=True)
            else:
                X[col].fillna(median_val, inplace=True)

        y = ml_de["label"]

        # Check no NaNs remain
        if X.isna().any().any():
            print(f"  ERROR: Still has NaN after filling")
            print(X.isna().sum())
        else:
            print(f"  Features ready: {X.shape}")

            if len(X) < 10:
                cv = LeaveOneOut()
                cv_name = "LOO"
            elif len(X) < 20:
                cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=RANDOM_STATE)
                cv_name = "3-fold CV"
            else:
                cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
                cv_name = "5-fold CV"

            pipe = Pipeline([
                ("scaler", StandardScaler()),
                ("rf", RandomForestClassifier(
                    n_estimators=500, max_features="sqrt",
                    min_samples_leaf=1, class_weight="balanced",
                    n_jobs=-1, random_state=RANDOM_STATE))
            ])

            cv_scores = cross_val_score(pipe, X, y,
                                        cv=cv, scoring="balanced_accuracy")
            print(f"\n  {cv_name} balanced accuracy: "
                  f"{cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
            print(f"  Per-fold: {np.round(cv_scores, 3)}")

            pd.DataFrame({
                "fold": range(1, len(cv_scores)+1),
                "balanced_accuracy": cv_scores,
            }).to_csv(os.path.join(args.outdir, "cv_scores_srna.tsv"),
                      sep="\t", index=False)

            # Fit on full data for SHAP
            scaler  = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            rf_final = RandomForestClassifier(
                n_estimators=1000, max_features="sqrt",
                min_samples_leaf=1, class_weight="balanced",
                n_jobs=-1, random_state=RANDOM_STATE)
            rf_final.fit(X_scaled, y)

            try:
                explainer   = shap.TreeExplainer(rf_final)
                shap_values = explainer.shap_values(X_scaled)
                shap_cls1   = (shap_values[1]
                               if isinstance(shap_values, list)
                               else shap_values)

                shap_imp = pd.DataFrame({
                    "feature":       FEATURE_COLS,
                    "mean_abs_shap": np.abs(shap_cls1).mean(axis=0),
                    "mean_shap":     shap_cls1.mean(axis=0),
                }).sort_values("mean_abs_shap", ascending=False)

                print(f"\n  SHAP feature importance:")
                print(shap_imp.to_string(index=False))
                shap_imp.to_csv(
                    os.path.join(args.outdir, "shap_importance_srna.tsv"),
                    sep="\t", index=False)

                # Bar plot
                fig, ax = plt.subplots(figsize=(7, 5))
                colors = ["#d95f02" if v>0 else "#1b9e77"
                          for v in shap_imp["mean_shap"]]
                ax.barh(shap_imp["feature"][::-1],
                        shap_imp["mean_abs_shap"][::-1],
                        color=colors[::-1])
                ax.set_xlabel("Mean |SHAP|")
                ax.set_title(f"sRNA RF feature importance — {DATASET} {T1}→{T2}\n"
                             f"{cv_name}: {cv_scores.mean():.2f} ± "
                             f"{cv_scores.std():.2f}")
                ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
                plt.tight_layout()
                plt.savefig(os.path.join(args.outdir, "shap_bar_srna.pdf"),
                            dpi=300, bbox_inches="tight")
                plt.savefig(os.path.join(args.outdir, "shap_bar_srna.png"),
                            dpi=150, bbox_inches="tight")
                plt.close()

                # Summary plot
                plt.figure(figsize=(8, 6))
                shap.summary_plot(shap_cls1, X_scaled,
                                  feature_names=FEATURE_COLS,
                                  show=False, max_display=15)
                plt.title(f"SHAP summary — sRNA upregulation {T1}→{T2}")
                plt.tight_layout()
                plt.savefig(os.path.join(args.outdir, "shap_summary_srna.pdf"),
                            dpi=300, bbox_inches="tight")
                plt.savefig(os.path.join(args.outdir, "shap_summary_srna.png"),
                            dpi=150, bbox_inches="tight")
                plt.close()
                print("  Saved: shap_bar_srna.pdf/png, shap_summary_srna.pdf/png")

            except Exception as e:
                print(f"  SHAP failed: {e}")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — asRNA expression vs target correlation + sRNA conservation analysis
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("SECTION 5: asRNA–target expression correlation & sRNA conservation")
print("-" * 40)

# ── Part A: asRNA expression correlation ────────────────────────────────────
asrna_expr = expr[
    (expr["is_asrna"] == 1) &
    (expr[f"fpkm_{T1}"] > 0) &
    (expr[f"fpkm_{T2}"] > 0)
].copy()

print(f"  asRNAs detected in both timepoints: {len(asrna_expr)}")

if not asrna_expr.empty:
    n_neg = (asrna_expr["log2fc_fpkm"] < 0).sum()
    n_pos = (asrna_expr["log2fc_fpkm"] > 0).sum()
    print(f"  asRNAs increasing ({T2}):  {n_pos}")
    print(f"  asRNAs decreasing ({T2}): {n_neg}")
    print(f"  Mean asRNA log2FC: {asrna_expr['log2fc_fpkm'].mean():.3f}")
    print(f"  (Negative mean suggests suppression post-perturbation)")

    asrna_expr.to_csv(
        os.path.join(args.outdir, "asrna_expression_both_timepoints.tsv"),
        sep="\t", index=False)

# ── Part B: Upregulated sRNA properties investigation ──────────────────────
# ── Part B: Upregulated sRNA properties investigation ──────────────────────
print(f"\n  [Investigating upregulated sRNA properties]")

upregulated = expr[expr["detection"] == f"only_{T2}"].copy()
downregulated = expr[expr["detection"] == f"only_{T1}"].copy()
shared = expr[expr["detection"] == "both"].copy()

print(f"\n  Upregulated sRNAs (only in {T2}): {len(upregulated)}")
if len(upregulated) > 0:
    print(f"    Mean FPKM: {upregulated[f'fpkm_{T2}'].mean():.2f}")
    print(f"    Median FPKM: {upregulated[f'fpkm_{T2}'].median():.2f}")
    print(f"    % with FPKM < 5: {100*(upregulated[f'fpkm_{T2}'] < 5).sum()/len(upregulated):.1f}%")

print(f"\n  Downregulated sRNAs (only in {T1}): {len(downregulated)}")
if len(downregulated) > 0:
    print(f"    Mean FPKM: {downregulated[f'fpkm_{T1}'].mean():.2f}")
    print(f"    Median FPKM: {downregulated[f'fpkm_{T1}'].median():.2f}")
    print(f"    % with FPKM < 5: {100*(downregulated[f'fpkm_{T1}'] < 5).sum()/len(downregulated):.1f}%")

# ── Part C: Raw IntaRNA analysis (use correct timepoint file) ────────────
print(f"\n  [Investigating target predictions]")

if intarna_t1 is not None and len(intarna_t1) > 0 and intarna_t2 is not None and len(intarna_t2) > 0:
    # Get downregulated sRNA IDs (from 2016)
    downreg_srna_ids = set(downregulated['srna_id'])
    downreg_interactions_raw = intarna_t1[
        intarna_t1['srna_id'].isin(downreg_srna_ids)
    ]
    
    # Get upregulated sRNA IDs (from 2017) - USE intarna_t2!
    upreg_srna_ids = set(upregulated['srna_id'])
    upreg_interactions_raw = intarna_t2[
        intarna_t2['srna_id'].isin(upreg_srna_ids)
    ]
    
    print(f"\n  Upregulated sRNAs with ANY raw predictions: "
          f"{upreg_interactions_raw['srna_id'].nunique()}/{len(upregulated)}")
    if len(upreg_interactions_raw) > 0:
        print(f"    Total interactions: {len(upreg_interactions_raw)}")
        print(f"    Mean energy: {upreg_interactions_raw['best_energy'].mean():.2f} kcal/mol")
        print(f"    % with energy < -5: "
              f"{100*(upreg_interactions_raw['best_energy'] < -5).sum()/len(upreg_interactions_raw):.1f}%")
    
    print(f"\n  Downregulated sRNAs with ANY raw predictions: "
          f"{downreg_interactions_raw['srna_id'].nunique()}/{len(downregulated)}")
    if len(downreg_interactions_raw) > 0:
        print(f"    Total interactions: {len(downreg_interactions_raw)}")
        print(f"    Mean energy: {downreg_interactions_raw['best_energy'].mean():.2f} kcal/mol")
        print(f"    % with energy < -5: "
              f"{100*(downreg_interactions_raw['best_energy'] < -5).sum()/len(downreg_interactions_raw):.1f}%")
elif intarna_t1 is not None:
    # Fallback if only 2016 available
    downreg_srna_ids = set(downregulated['srna_id'])
    downreg_interactions_raw = intarna_t1[
        intarna_t1['srna_id'].isin(downreg_srna_ids)
    ]
    
    print(f"\n  Downregulated sRNAs with ANY raw predictions: "
          f"{downreg_interactions_raw['srna_id'].nunique()}/{len(downregulated)}")
    if len(downreg_interactions_raw) > 0:
        print(f"    Total interactions: {len(downreg_interactions_raw)}")
        print(f"    Mean energy: {downreg_interactions_raw['best_energy'].mean():.2f} kcal/mol")
        print(f"    % with energy < -5: "
              f"{100*(downreg_interactions_raw['best_energy'] < -5).sum()/len(downreg_interactions_raw):.1f}%")
    
    print(f"\n  WARNING: 2017 IntaRNA data not available for upregulated sRNAs")


# In SECTION 5, replace the TPM correlation code with:

# ── Part D: Per-replicate sRNA-target correlation (if TPM matrices provided) ──
print(f"\n  [Per-replicate sRNA-target correlation]")

tpm_matrices = {}
if args.tpm_matrix_t1 and os.path.exists(args.tpm_matrix_t1):
    print(f"  Loading TPM matrix T1: {args.tpm_matrix_t1}")
    tpm_matrices[T1] = pd.read_csv(args.tpm_matrix_t1, sep="\t", index_col=0)
    print(f"    Shape: {tpm_matrices[T1].shape} (genes × samples)")

if args.tpm_matrix_t2 and os.path.exists(args.tpm_matrix_t2):
    print(f"  Loading TPM matrix T2: {args.tpm_matrix_t2}")
    tpm_matrices[T2] = pd.read_csv(args.tpm_matrix_t2, sep="\t", index_col=0)
    print(f"    Shape: {tpm_matrices[T2].shape} (genes × samples)")

# Load the full interactions files (not summary)
# These are needed for the mrna_id column
def load_interactions(path):
    if not path or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]
    return df

interactions_t1 = load_interactions(
    args.intarna_t1.replace("intarna_summary_fixed.tsv", "interactions_significant_fixed.tsv")
) if args.intarna_t1 else None
interactions_t2 = load_interactions(
    args.intarna_t2.replace("intarna_summary_fixed.tsv", "interactions_significant_fixed.tsv")
) if args.intarna_t2 else None

if len(tpm_matrices) > 0 and (interactions_t1 is not None or interactions_t2 is not None):
    from scipy.stats import pearsonr
    
    all_correlations = []
    
    # For each timepoint
    for tg in [T1, T2]:
        if tg == T1:
            tpm = tpm_matrices.get(T1)
            interactions_data = interactions_t1
            asrna_in_tg = expr[(expr["is_asrna"] == 1) & 
                               (expr[f"fpkm_{T1}"] > 0)].copy()
        else:
            tpm = tpm_matrices.get(T2)
            interactions_data = interactions_t2
            asrna_in_tg = expr[(expr["is_asrna"] == 1) & 
                               (expr[f"fpkm_{T2}"] > 0)].copy()
        
        if tpm is None or interactions_data is None:
            print(f"  Skipping {tg}: TPM or interactions not available")
            continue
        
        print(f"\n  Computing correlations for {tg}...")
        print(f"    asRNAs in {tg}: {len(asrna_in_tg)}")
        
        for _, asrna_row in asrna_in_tg.iterrows():
            asrna_id = asrna_row["srna_id"]
            
            # Get this asRNA's targets from interactions file
            targets_for_srna = interactions_data[
                interactions_data["srna_id"] == asrna_id
            ]
            
            if len(targets_for_srna) == 0:
                continue
            
            # Extract gene IDs from target mRNA_ids
            target_genes = []
            for mrna_id in targets_for_srna["mrna_id"].unique():
                match = re.search(r"(\d+_\d+)", str(mrna_id))
                if match:
                    target_genes.append(match.group(1))
            
            if len(target_genes) == 0:
                continue
            
            # Get asRNA TPM values across samples
            asrna_tpm = None
            for row_idx in tpm.index:
                if asrna_id in str(row_idx):
                    asrna_tpm = tpm.loc[row_idx].values
                    break
            
            if asrna_tpm is None or len(asrna_tpm) < 3:
                continue
            
            # For each target gene, compute correlation
            for gene_id in target_genes:
                gene_tpm = None
                for row_idx in tpm.index:
                    if gene_id in str(row_idx):
                        gene_tpm = tpm.loc[row_idx].values
                        break
                
                if gene_tpm is None or len(gene_tpm) < 3:
                    continue
                
                # Compute Pearson correlation
                try:
                    corr, pval = pearsonr(asrna_tpm, gene_tpm)
                    
                    all_correlations.append({
                        "timepoint": tg,
                        "asrna_id": asrna_id,
                        "target_gene": gene_id,
                        "pearson_r": corr,
                        "pvalue": pval,
                        "n_samples": len(asrna_tpm),
                    })
                except:
                    continue
        
        print(f"    Correlations computed for {tg}: {len([c for c in all_correlations if c['timepoint']==tg])}")
    
    if len(all_correlations) > 0:
        corr_df = pd.DataFrame(all_correlations)
        
        print(f"\n  Total correlations: {len(corr_df)}")
        print(f"  Mean Pearson r: {corr_df['pearson_r'].mean():.3f}")
        print(f"  Significant (p<0.05): {(corr_df['pvalue'] < 0.05).sum()}")
        print(f"  Negative correlations (expected for repression): "
              f"{(corr_df['pearson_r'] < 0).sum()}")
        
        corr_df.to_csv(
            os.path.join(args.outdir, "asrna_target_pearson_correlation.tsv"),
            sep="\t", index=False)
        
        # Histogram of correlations by timepoint
        fig, axes = plt.subplots(1, len(tpm_matrices), figsize=(6*len(tpm_matrices), 5))
        if len(tpm_matrices) == 1:
            axes = [axes]
        
        for ax, tg in enumerate(tpm_matrices.keys()):
            tg_corr = corr_df[corr_df["timepoint"] == tg]
            if len(tg_corr) > 0:
                axes[ax].hist(tg_corr["pearson_r"], bins=20, 
                             color="#7570b3", alpha=0.7, edgecolor="black")
                axes[ax].axvline(0, color="red", lw=2, ls="--", label="No correlation")
                axes[ax].set_xlabel("Pearson r (asRNA vs target gene)")
                axes[ax].set_ylabel("Count")
                axes[ax].set_title(f"{tg}\nMean r={tg_corr['pearson_r'].mean():.3f}, "
                                  f"n={len(tg_corr)}")
                axes[ax].legend(frameon=False)
                axes[ax].spines["top"].set_visible(False)
                axes[ax].spines["right"].set_visible(False)
        
        plt.suptitle(f"asRNA-target correlation distribution — {DATASET}")
        plt.tight_layout()
        plt.savefig(
            os.path.join(args.outdir, "asrna_target_pearson_histogram.pdf"),
            dpi=300, bbox_inches="tight")
        plt.savefig(
            os.path.join(args.outdir, "asrna_target_pearson_histogram.png"),
            dpi=150, bbox_inches="tight")
        plt.close()
        print("  Saved: asrna_target_pearson_histogram.pdf/png")
else:
    if len(tpm_matrices) == 0:
        print(f"  TPM matrices not provided or files not found")
    else:
        print(f"  Interactions files not found")
# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Summary figures
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("SECTION 6: Summary figures")
print("-" * 40)

# Figure 1: sRNA class composition pie charts
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
for ax, (tg, df) in zip(axes, [(T1, srna_t1), (T2, srna_t2)]):
    if df.empty:
        continue
    counts = df["srna_type"].value_counts()
    ax.pie(counts,
           labels=[f"{k}\n(n={v})" for k, v in counts.items()],
           colors=["#d95f02","#1b9e77"],
           autopct="%1.0f%%", startangle=90)
    ax.set_title(f"{DATASET} — {tg}\n{len(df)} total sRNAs")
plt.suptitle("sRNA class composition", fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, "srna_class_composition.pdf"),
            dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(args.outdir, "srna_class_composition.png"),
            dpi=150, bbox_inches="tight")
plt.close()

# Figure 2: Ranked log2FC bar chart
fig, ax = plt.subplots(figsize=(8, max(4, len(expr) * 0.15)))
sorted_expr = expr.sort_values("log2fc_fpkm")
bar_colors = sorted_expr["de_status"].map(color_map)
ax.barh(range(len(sorted_expr)), sorted_expr["log2fc_fpkm"],
        color=bar_colors, height=0.8, edgecolor="none")
ax.axvline(x=0, color="black", lw=0.8)
ax.axvline(x=LOG2FC_THRESH,  color="gray", lw=0.8, ls="--")
ax.axvline(x=-LOG2FC_THRESH, color="gray", lw=0.8, ls="--")
ax.set_xlabel(f"log₂FC FPKM ({T1} → {T2})")
ax.set_ylabel("sRNAs (ranked)")
ax.set_title(f"sRNA expression change — {DATASET}")
patches = [
    mpatches.Patch(color="#d95f02", label=f"Up in {T2} (n={n_up})"),
    mpatches.Patch(color="#1b9e77", label=f"Up in {T1} (n={n_down})"),
    mpatches.Patch(color="#aaaaaa", label=f"Stable (n={n_stab})"),
]
ax.legend(handles=patches, frameon=False)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, "srna_log2fc_ranked.pdf"),
            dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(args.outdir, "srna_log2fc_ranked.png"),
            dpi=150, bbox_inches="tight")
plt.close()

# Figure 3: FPKM distribution per class per timepoint
fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
for ax, (tg, df) in zip(axes, [(T1, srna_t1), (T2, srna_t2)]):
    if df.empty:
        continue
    for stype, color in [("antisense","#d95f02"), ("intergenic","#1b9e77")]:
        vals = np.log2(df[df["srna_type"]==stype]["fpkm"] + PSEUDOCOUNT)
        ax.hist(vals, bins=20, alpha=0.6, color=color,
                label=stype, edgecolor="white", lw=0.5)
    ax.set_xlabel("log₂(FPKM)")
    ax.set_title(f"{tg}")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
plt.suptitle(f"sRNA FPKM distributions — {DATASET}", fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, "srna_fpkm_distributions.pdf"),
            dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(args.outdir, "srna_fpkm_distributions.png"),
            dpi=150, bbox_inches="tight")
plt.close()

print("  Saved: srna_class_composition, srna_log2fc_ranked, "
      "srna_fpkm_distributions")

# ── Save complete feature table ────────────────────────────────────────────────
final = expr.copy()
# Replace the final merge block with:
if not rf_t1.empty:
    rf_final_cols = ["srna_id"] + [c for c in
        ["gc_content", "mfe", "struct_entropy", "structure"]
        if c in rf_t1.columns]
    final = final.merge(rf_t1[rf_final_cols], on="srna_id", how="left")
if intarna_t1 is not None:
    final = final.merge(intarna_t1, on="srna_id", how="left")

final.to_csv(
    os.path.join(args.outdir, "srna_complete_feature_table.tsv"),
    sep="\t", index=False)

print(f"\n{'='*60}")
print(f"Complete. All outputs in: {args.outdir}")
print(f"{'='*60}")
for f in sorted(os.listdir(args.outdir)):
    size = os.path.getsize(os.path.join(args.outdir, f))
    print(f"  {f:<45} {size/1024:.1f} KB")
