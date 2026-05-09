#!/usr/bin/env python3
"""
07_build_feature_matrix.py

Builds the Random Forest input feature matrix.
Primary expression feature: FPKM log2FC between timegroups from SnapT GFF.
IntaRNA is optional — included if available, skipped gracefully if not.

Usage:
    python 07_build_feature_matrix.py \
        --dataset    AD \
        --metadata   /path/to/Metadata.csv \
        --snapt      /path/to/AD/snapt_2016/small_antisense_ncRNAs.gff \
        --rnafold    /path/to/AD/rnafold_2016/rnafold_summary.tsv \
        --humann_rna /path/to/AD/Total_RNA_pathabundance_humann_table.tsv \
        --humann_dna /path/to/AD/Total_DNA_pathabundance_humann_table.tsv \
        --metaphlan  /path/to/AD/Total_DNA_merged_abundance_table.tsv \
        --output     /path/to/AD/feature_matrix.tsv \
        [--intarna   /path/to/AD/intarna_2016/intarna_summary.tsv]
"""

import argparse
import json
import numpy as np
import pandas as pd
import re
import os

# ── Arguments ─────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--dataset",    required=True,  help="AD or CB")
parser.add_argument("--metadata",   required=True)
parser.add_argument("--snapt",      required=True,  help="Path to any SnapT GFF — used to locate dataset dir")
parser.add_argument("--rnafold",    required=True,  help="RNAfold summary TSV for first timegroup")
parser.add_argument("--humann_rna", required=True)
parser.add_argument("--humann_dna", required=True)
parser.add_argument("--metaphlan",  required=True)
parser.add_argument("--output",     required=True)
parser.add_argument("--intarna",    required=False, default=None,
                    help="IntaRNA summary TSV (optional)")
args = parser.parse_args()

PSEUDOCOUNT = 1e-6
DATASET     = args.dataset.upper()
assert DATASET in ("AD", "CB"), "Dataset must be AD or CB"

# Dataset base directory — parent of the snapt_* directories
DATASET_DIR = os.path.dirname(os.path.dirname(args.snapt))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Metadata and timepoints
# ═══════════════════════════════════════════════════════════════════════════════
print(f"Loading metadata: {DATASET}")

meta = pd.read_csv(args.metadata, sep="\t")
print(f"  Metadata columns: {list(meta.columns)}")
print(f"  First row: {meta.iloc[0].to_dict()}")
# ─────────────────────────────────────────────────────────────────────────────
meta.columns = [c.strip() for c in meta.columns]
meta = meta[meta["Dataset"].str.strip() == DATASET].copy()

rna_meta = meta[meta["Type"] == "Metatranscriptome"][["Name"]].rename(columns={"Name": "sample_id"})
dna_meta = meta[meta["Type"] == "Metagenome"][["Name"]].rename(columns={"Name": "sample_id"})

print(f"  RNA samples: {len(rna_meta)}  DNA samples: {len(dna_meta)}")

def extract_timepoint(name):
    if DATASET == "AD":
        m = re.search(r"_(201[0-9])_", str(name))
        return m.group(1) if m else "unknown"
    else:
        m = re.search(r"CB_(Fall|Sum)_", str(name))
        return m.group(1) if m else "unknown"

rna_meta["timepoint"] = rna_meta["sample_id"].apply(extract_timepoint)
dna_meta["timepoint"] = dna_meta["sample_id"].apply(extract_timepoint)

timepoints = sorted(rna_meta["timepoint"].unique())
tp_pairs   = [(timepoints[i], timepoints[i+1]) for i in range(len(timepoints)-1)]
print(f"  Timepoints: {timepoints}")
print(f"  Pairs:      {tp_pairs}")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Parse SnapT GFF files for all timegroups
#
# Key insight: FPKM and TPM in the SnapT GFF are direct expression estimates
# from StringTie — more reliable than derived counts and available immediately.
# log2FC of FPKM between timegroups is our primary expression feature.
# ═══════════════════════════════════════════════════════════════════════════════
print("\nParsing SnapT GFF files...")

def parse_snapt_gff(gff_path, srna_type):
    """Parse SnapT GFF into DataFrame. Returns empty DF if file missing."""
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

            transcript_id = get_attr("transcript_id")
            cov           = float(get_attr("cov")  or 0)
            fpkm          = float(get_attr("FPKM") or 0)
            tpm           = float(get_attr("TPM")  or 0)
            antisense_to  = get_attr("antisense_to_gene")
            srna_id       = f"{contig}:{start-1}-{end}({strand})"

            rows.append({
                "srna_id":       srna_id,
                "transcript_id": transcript_id,
                "contig":        contig,
                "start":         start,
                "end":           end,
                "strand":        strand,
                "length_gff":    end - start + 1,
                "srna_type":     srna_type,
                "is_asrna":      1 if srna_type == "antisense" else 0,
                "cov":           cov,
                "fpkm":          fpkm,
                "tpm":           tpm,
                "antisense_to":  antisense_to,
            })
    return pd.DataFrame(rows)

# Load GFFs for each timegroup
snapt_by_tg = {}
for tg in timepoints:
    tg_dir     = os.path.join(DATASET_DIR, f"snapt_{tg}")
    asrna_gff  = os.path.join(tg_dir, "small_antisense_ncRNAs.gff")
    itsrna_gff = os.path.join(tg_dir, "small_intergenic_ncRNAs.gff")

    frames = []
    for path, stype in [(asrna_gff, "antisense"), (itsrna_gff, "intergenic")]:
        df = parse_snapt_gff(path, stype)
        if not df.empty:
            df["timegroup"] = tg
            frames.append(df)
            print(f"  {tg} {stype}: {len(df)} sRNAs  "
                  f"(FPKM mean={df['fpkm'].mean():.1f})")

    snapt_by_tg[tg] = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

# All sRNAs combined
snapt_all = pd.concat(snapt_by_tg.values(), ignore_index=True)
print(f"  Total sRNAs: {len(snapt_all)}  "
      f"(asRNA={snapt_all['is_asrna'].sum()}  "
      f"itsRNA={(snapt_all['is_asrna']==0).sum()})")

# ── Compute FPKM and TPM log2FC between consecutive timegroups ────────────────
# This is the PRIMARY expression feature replacing IntaRNA-derived features.
# For each sRNA present in both timegroups, compute log2FC.
# sRNAs only in one timegroup get log2FC = log2(fpkm+pseudo / pseudo) or vice versa.

print("\nComputing FPKM log2FC between timegroups...")
log2fc_frames = []

for t1, t2 in tp_pairs:
    df1 = snapt_by_tg.get(t1, pd.DataFrame())
    df2 = snapt_by_tg.get(t2, pd.DataFrame())

    if df1.empty or df2.empty:
        print(f"  WARNING: Missing SnapT data for {t1} or {t2} — skipping log2FC")
        continue

    # Pivot to get FPKM per srna_id per timegroup
    fpkm1 = df1.set_index("srna_id")["fpkm"].rename(f"fpkm_{t1}")
    fpkm2 = df2.set_index("srna_id")["fpkm"].rename(f"fpkm_{t2}")
    tpm1  = df1.set_index("srna_id")["tpm"].rename(f"tpm_{t1}")
    tpm2  = df2.set_index("srna_id")["tpm"].rename(f"tpm_{t2}")

    # Union of all sRNA IDs across both timegroups
    all_ids = fpkm1.index.union(fpkm2.index)
    merged  = pd.DataFrame(index=all_ids)
    merged  = merged.join(fpkm1).join(fpkm2).join(tpm1).join(tpm2).fillna(0)

    merged[f"log2fc_fpkm_{t1}_vs_{t2}"] = np.log2(
        (merged[f"fpkm_{t2}"] + PSEUDOCOUNT) /
        (merged[f"fpkm_{t1}"] + PSEUDOCOUNT)
    )
    merged[f"log2fc_tpm_{t1}_vs_{t2}"] = np.log2(
        (merged[f"tpm_{t2}"]  + PSEUDOCOUNT) /
        (merged[f"tpm_{t1}"]  + PSEUDOCOUNT)
    )
    merged["pair_label"] = f"{t1}_vs_{t2}"
    merged = merged.reset_index().rename(columns={"index": "srna_id"})
    log2fc_frames.append(merged)
    print(f"  {t1}→{t2}: {len(merged)} sRNAs  "
          f"log2FC FPKM range "
          f"[{merged[f'log2fc_fpkm_{t1}_vs_{t2}'].min():.2f}, "
          f"{merged[f'log2fc_fpkm_{t1}_vs_{t2}'].max():.2f}]")

srna_log2fc = pd.concat(log2fc_frames, ignore_index=True) if log2fc_frames else pd.DataFrame()

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — RNAfold structural features
# ═══════════════════════════════════════════════════════════════════════════════
print("\nLoading RNAfold results...")

rnafold_frames = []
for tg in timepoints:
    rf_path = os.path.join(DATASET_DIR, f"rnafold_{tg}", "rnafold_summary.tsv")
    if os.path.exists(rf_path):
        df = pd.read_csv(rf_path, sep="\t")
        df["timegroup"] = tg
        rnafold_frames.append(df)
        print(f"  {tg}: {len(df)} structures loaded")
    else:
        print(f"  WARNING: {rf_path} not found — RNAfold features will be NaN for {tg}")

rnafold = pd.concat(rnafold_frames, ignore_index=True) if rnafold_frames else pd.DataFrame()

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — IntaRNA (timepoint-aware, like RNAfold)
# ═══════════════════════════════════════════════════════════════════════════════
print("\nLoading IntaRNA results...")

intarna_frames = []

for tg in timepoints:
    ia_path = os.path.join(DATASET_DIR, f"intarna_{tg}", "intarna_summary.tsv")

    if os.path.exists(ia_path):
        df = pd.read_csv(ia_path, sep="\t")
        df.columns = [c.strip().lower() for c in df.columns]
        df["timegroup"] = tg
        intarna_frames.append(df)

        print(f"  {tg}: {len(df)} sRNAs "
              f"(mean energy={df['best_energy'].mean():.2f})")
    else:
        print(f"  WARNING: {ia_path} not found — IntaRNA features NaN for {tg}")

intarna = (
    pd.concat(intarna_frames, ignore_index=True)
    if intarna_frames else pd.DataFrame()
)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Merge all sRNA features per sRNA
# ═══════════════════════════════════════════════════════════════════════════════
print("\nMerging sRNA feature tables...")

# Base: all sRNAs with structural features
srna = snapt_all.copy()

if not rnafold.empty and "srna_id" in rnafold.columns:
    srna = srna.merge(
        rnafold[["srna_id", "gc_content", "mfe", "ensemble_energy", "struct_entropy"]],
        on="srna_id", how="left"
    )

if intarna is not None and "srna_id" in intarna.columns:
    srna = srna.merge(
        intarna,
        on=["srna_id", "timegroup"],
        how="left"
    )

    for col in ["n_targets", "best_energy", "mean_energy", "std_energy"]:
        if col in srna.columns:
            srna[col] = srna[col].fillna(0)

print(f"  sRNA feature table: {srna.shape}")

print("\nComputing IntaRNA feature changes...")

intarna_delta_frames = []

for t1, t2 in tp_pairs:
    df1 = intarna[intarna["timegroup"] == t1]
    df2 = intarna[intarna["timegroup"] == t2]

    if df1.empty or df2.empty:
        continue

    f1 = df1.set_index("srna_id")
    f2 = df2.set_index("srna_id")

    all_ids = f1.index.union(f2.index)
    merged = pd.DataFrame(index=all_ids)

    for col in ["n_targets", "best_energy", "mean_energy"]:
        if col in f1.columns and col in f2.columns:
            merged[f"{col}_{t1}"] = f1[col]
            merged[f"{col}_{t2}"] = f2[col]

            merged[f"delta_{col}_{t1}_vs_{t2}"] = (
                merged[f"{col}_{t2}"].fillna(0) -
                merged[f"{col}_{t1}"].fillna(0)
            )

    merged["pair_label"] = f"{t1}_vs_{t2}"
    merged = merged.reset_index().rename(columns={"index": "srna_id"})
    intarna_delta_frames.append(merged)

intarna_delta = (
    pd.concat(intarna_delta_frames, ignore_index=True)
    if intarna_delta_frames else pd.DataFrame()
)
# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — HUMAnN expression proxy
# ═══════════════════════════════════════════════════════════════════════════════
print("\nLoading HUMAnN outputs...")

def normalize_col(name):
    return re.sub(
        r"(_paired_combined|_paired|_combined"
        r"|_knead_out|_kneaddata"
        r"|_humann|_pathabundance|_Abundance"
        r"|_profile|_RNA\d*|_DNA).*",
        "", name, flags=re.IGNORECASE
    )

def load_humann(path):
    # Read skipping comment lines but keeping the actual header
    df = pd.read_csv(path, sep="\t", comment="#", header=0)
    # If comment="#" ate the header, re-read without it
    if df.columns[0].startswith("#"):
        df.columns = [c.lstrip("# ") for c in df.columns]
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={df.columns[0]: "Pathway"})
    df = df[~df["Pathway"].str.contains("UNMAPPED|UNINTEGRATED", na=False)]
    df = df[~df["Pathway"].str.contains(r"\|", na=False)]
    df.columns = ["Pathway"] + [normalize_col(c) for c in df.columns[1:]]
    return df

humann_rna = load_humann(args.humann_rna)
humann_dna = load_humann(args.humann_dna)

rna_activity = humann_rna.set_index("Pathway").sum(axis=0).rename("rna_activity")
dna_activity = humann_dna.set_index("Pathway").sum(axis=0).rename("dna_activity")

activity_df = pd.DataFrame({
    "rna_activity": rna_activity,
    "dna_activity": dna_activity,
}).reset_index().rename(columns={"index": "sample_id"})

activity_df["sample_id_norm"] = activity_df["sample_id"].apply(normalize_col)
rna_meta["sample_id_norm"]    = rna_meta["sample_id"].apply(normalize_col)

activity_df = activity_df.merge(
    rna_meta[["sample_id_norm", "timepoint"]], on="sample_id_norm", how="left"
)

tp_activity = (
    activity_df.groupby("timepoint")
    .agg(
        mean_rna_activity = ("rna_activity", "mean"),
        std_rna_activity  = ("rna_activity", "std"),
        mean_dna_activity = ("dna_activity", "mean"),
    )
    .reset_index()
)
print(f"  Timepoint activity:\n{tp_activity.to_string(index=False)}")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — MetaPhlAn abundance (target variable)
# ═══════════════════════════════════════════════════════════════════════════════
print("\nLoading MetaPhlAn abundance...")

mpa = pd.read_csv(args.metaphlan, sep="\t", comment="#")
mpa.columns = [c.strip() for c in mpa.columns]
mpa = mpa.rename(columns={mpa.columns[0]: "clade_name"})
mpa = mpa[
    mpa["clade_name"].str.contains(r"s__", na=False) &
    ~mpa["clade_name"].str.contains(r"t__", na=False)
].copy()

# ── FIX: rename sample columns first, then add species ───────────────────────
sample_cols = mpa.columns[1:]  # everything except clade_name
mpa.columns = ["clade_name"] + [normalize_col(c) for c in sample_cols]
mpa["species"] = mpa["clade_name"].str.extract(r"s__([^\|]+)$")
# ─────────────────────────────────────────────────────────────────────────────

mpa_long = mpa.melt(
    id_vars=["clade_name", "species"],
    var_name="sample_id_norm", value_name="rel_abundance"
)

all_meta = pd.concat([rna_meta, dna_meta])
all_meta["sample_id_norm"] = all_meta["sample_id"].apply(normalize_col)

mpa_long = mpa_long.merge(
    all_meta[["sample_id_norm", "timepoint"]].drop_duplicates(),
    on="sample_id_norm", how="left"
)
mpa_long["rel_abundance"] = pd.to_numeric(mpa_long["rel_abundance"], errors="coerce")

tp_abundance = (
    mpa_long.dropna(subset=["timepoint"])
    .groupby(["species", "timepoint"])
    .agg(mean_abundance=("rel_abundance", "mean"))
    .reset_index()
)
print(f"  Species detected: {tp_abundance['species'].nunique()}")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — Build feature matrix
# ═══════════════════════════════════════════════════════════════════════════════
print("\nBuilding feature matrix...")

# Columns to aggregate per timepoint pair
SRNA_STATIC = [c for c in [
    "length_gff", "gc_content", "mfe", "ensemble_energy",
    "struct_entropy", "cov", "fpkm", "tpm", "is_asrna",
] if c in srna.columns]

INTARNA_COLS = [c for c in [
    "n_targets", "best_energy", "mean_energy", "std_energy"
] if c in srna.columns]

timepoint_rows = []
species_rows   = []

for t1, t2 in tp_pairs:
    print(f"  Pair: {t1} → {t2}")
    pair_label = f"{t1}_vs_{t2}"

    agg = {
        "dataset":      DATASET,
        "timepoint_T1": t1,
        "timepoint_T2": t2,
        "pair_label":   pair_label,
    }

    # ── Static sRNA structural/sequence features ──────────────────────────────
    # Use T1 timegroup sRNAs as the baseline regulatory state
    srna_t1 = srna[srna["timegroup"] == t1] if "timegroup" in srna.columns else srna

    n_total = len(srna_t1)
    n_asrna = srna_t1["is_asrna"].sum() if n_total else 0

    for col in SRNA_STATIC:
        vals = srna_t1[col].dropna()
        agg[f"srna_mean_{col}"] = vals.mean() if len(vals) else np.nan
        agg[f"srna_std_{col}"]  = vals.std()  if len(vals) else np.nan

    agg["asrna_fraction"]     = n_asrna / n_total if n_total else np.nan
    agg["itsrna_fraction"]    = 1 - agg["asrna_fraction"]
    agg["asrna_itsrna_ratio"] = n_asrna / (n_total - n_asrna + PSEUDOCOUNT)
    agg["n_srnas_T1"]         = n_total

    # ── IntaRNA features (if available) ───────────────────────────────────────
    for col in INTARNA_COLS:
        vals = srna_t1[col].dropna()
        agg[f"srna_mean_{col}"] = vals.mean() if len(vals) else np.nan

    # ── FPKM log2FC — PRIMARY EXPRESSION FEATURE ─────────────────────────────
    # Aggregate across all sRNAs: how much did the sRNA pool change in expression?
    log2fc_col = f"log2fc_fpkm_{t1}_vs_{t2}"
    if not srna_log2fc.empty and log2fc_col in srna_log2fc.columns:
        pair_fc = srna_log2fc[srna_log2fc["pair_label"] == pair_label]
        fc_vals = pair_fc[log2fc_col].dropna()

        agg["mean_srna_log2fc_fpkm"]     = fc_vals.mean()
        agg["std_srna_log2fc_fpkm"]      = fc_vals.std()
        agg["n_srnas_upregulated"]        = (fc_vals > 0.5).sum()
        agg["n_srnas_downregulated"]      = (fc_vals < -0.5).sum()
        agg["asrna_mean_log2fc_fpkm"]    = (
            pair_fc.loc[pair_fc["srna_id"].isin(
                srna_t1[srna_t1["is_asrna"]==1]["srna_id"]
            ), log2fc_col].mean()
        )
        agg["itsrna_mean_log2fc_fpkm"]   = (
            pair_fc.loc[pair_fc["srna_id"].isin(
                srna_t1[srna_t1["is_asrna"]==0]["srna_id"]
            ), log2fc_col].mean()
        )
    else:
        for k in ["mean_srna_log2fc_fpkm", "std_srna_log2fc_fpkm",
                  "n_srnas_upregulated", "n_srnas_downregulated",
                  "asrna_mean_log2fc_fpkm", "itsrna_mean_log2fc_fpkm"]:
            agg[k] = np.nan

    # ── HUMAnN RNA:DNA ratio ───────────────────────────────────────────────────
    t1_act = tp_activity[tp_activity["timepoint"] == t1]
    t2_act = tp_activity[tp_activity["timepoint"] == t2]

    rna_t1 = t1_act["mean_rna_activity"].values[0] if len(t1_act) else np.nan
    dna_t1 = t1_act["mean_dna_activity"].values[0] if len(t1_act) else np.nan
    rna_t2 = t2_act["mean_rna_activity"].values[0] if len(t2_act) else PSEUDOCOUNT

    agg["rna_activity_T1"]     = rna_t1
    agg["rna_activity_std_T1"] = t1_act["std_rna_activity"].values[0] if len(t1_act) else np.nan
    agg["dna_activity_T1"]     = dna_t1
    agg["rna_dna_ratio_T1"]    = rna_t1 / (dna_t1 + PSEUDOCOUNT)
    agg["log2fc_rna_activity"] = np.log2(
        (rna_t2 + PSEUDOCOUNT) / (rna_t1 + PSEUDOCOUNT)
    )

    # ── Target variable: species abundance log2FC ─────────────────────────────
    t1_ab  = tp_abundance[tp_abundance["timepoint"]==t1].set_index("species")["mean_abundance"]
    t2_ab  = tp_abundance[tp_abundance["timepoint"]==t2].set_index("species")["mean_abundance"]
    shared = t1_ab.index.intersection(t2_ab.index)

    if len(shared):
        log2fc = np.log2(
            (t2_ab[shared].values + PSEUDOCOUNT) /
            (t1_ab[shared].values + PSEUDOCOUNT)
        )
        agg["target_mean_log2fc"]   = log2fc.mean()
        agg["target_std_log2fc"]    = log2fc.std()
        agg["n_species_increasing"] = (log2fc > 0).sum()
        agg["n_species_decreasing"] = (log2fc < 0).sum()
        agg["n_species_total"]      = len(shared)

        # Residual log2FC — regress out baseline abundance bias
        baseline = np.log10(t1_ab[shared].values + PSEUDOCOUNT)
        if np.std(baseline) > 0:
            coef      = np.polyfit(baseline, log2fc, 1)
            residuals = log2fc - np.polyval(coef, baseline)
            agg["target_residual_mean_log2fc"] = residuals.mean()
            agg["target_residual_std_log2fc"]  = residuals.std()
        else:
            agg["target_residual_mean_log2fc"] = agg["target_mean_log2fc"]
            agg["target_residual_std_log2fc"]  = agg["target_std_log2fc"]
    else:
        for k in ["target_mean_log2fc", "target_std_log2fc",
                  "n_species_increasing", "n_species_decreasing",
                  "n_species_total", "target_residual_mean_log2fc",
                  "target_residual_std_log2fc"]:
            agg[k] = np.nan

    timepoint_rows.append(agg)

    # ── IntaRNA delta features (interaction rewiring) ─────────────────────────
    delta_energy_col = f"delta_best_energy_{t1}_vs_{t2}"
    delta_targets_col = f"delta_n_targets_{t1}_vs_{t2}"

    if not intarna_delta.empty:
        pair_delta = intarna_delta[
            intarna_delta["pair_label"] == pair_label
        ]

        # Restrict to sRNAs present in T1 baseline (same logic as other features)
        pair_delta = pair_delta[
        pair_delta["srna_id"].isin(srna_t1["srna_id"])
        ]

        # ── Energy change ──────────────────────────────────────────────────────
        if delta_energy_col in pair_delta.columns:
            vals = pair_delta[delta_energy_col].dropna()
            agg["mean_delta_best_energy"] = vals.mean() if len(vals) else np.nan
            agg["std_delta_best_energy"]  = vals.std()  if len(vals) else np.nan
        else:
            agg["mean_delta_best_energy"] = np.nan
            agg["std_delta_best_energy"]  = np.nan

        # ── Target count change ────────────────────────────────────────────────
        if delta_targets_col in pair_delta.columns:
            vals = pair_delta[delta_targets_col].dropna()
            agg["mean_delta_n_targets"] = vals.mean() if len(vals) else np.nan
            agg["std_delta_n_targets"]  = vals.std()  if len(vals) else np.nan
        else:
            agg["mean_delta_n_targets"] = np.nan
            agg["std_delta_n_targets"]  = np.nan

        # ── Optional: directional summaries (VERY useful biologically) ─────────
        if delta_targets_col in pair_delta.columns:
            vals = pair_delta[delta_targets_col].fillna(0)
            agg["n_srnas_gained_targets"] = (vals > 0).sum()
            agg["n_srnas_lost_targets"]   = (vals < 0).sum()

    else:
        # No IntaRNA at all
        for k in [
            "mean_delta_best_energy", "std_delta_best_energy",
            "mean_delta_n_targets", "std_delta_n_targets",
            "n_srnas_gained_targets", "n_srnas_lost_targets"
        ]:
            agg[k] = np.nan
        
    # ── Species-level rows ────────────────────────────────────────────────────
    for sp in shared:
        species_rows.append({
            "dataset":        DATASET,
            "pair_label":     pair_label,
            "timepoint_T1":   t1,
            "timepoint_T2":   t2,
            "species":        sp,
            "abundance_T1":   float(t1_ab[sp]),
            "abundance_T2":   float(t2_ab[sp]),
            "log2fc":         float(np.log2((t2_ab[sp]+PSEUDOCOUNT)/(t1_ab[sp]+PSEUDOCOUNT))),
            "log10_baseline": float(np.log10(t1_ab[sp]+PSEUDOCOUNT)),
            # Attach timepoint-level sRNA features for species-level RF
            "mean_srna_log2fc_fpkm":   agg.get("mean_srna_log2fc_fpkm", np.nan),
            "asrna_fraction":          agg.get("asrna_fraction", np.nan),
            "asrna_itsrna_ratio":      agg.get("asrna_itsrna_ratio", np.nan),
            "rna_dna_ratio_T1":        agg.get("rna_dna_ratio_T1", np.nan),
            "log2fc_rna_activity":     agg.get("log2fc_rna_activity", np.nan),
            "n_srnas_upregulated":     agg.get("n_srnas_upregulated", np.nan),
            "n_srnas_downregulated":   agg.get("n_srnas_downregulated", np.nan),
            "asrna_mean_log2fc_fpkm":  agg.get("asrna_mean_log2fc_fpkm", np.nan),
            "itsrna_mean_log2fc_fpkm": agg.get("itsrna_mean_log2fc_fpkm", np.nan),
            "srna_mean_mfe":           agg.get("srna_mean_mfe", np.nan),
            "srna_mean_fpkm":          agg.get("srna_mean_fpkm", np.nan),
        })

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Save outputs
# ═══════════════════════════════════════════════════════════════════════════════
feat_matrix    = pd.DataFrame(timepoint_rows)
species_matrix = pd.DataFrame(species_rows)

feat_matrix.to_csv(args.output, sep="\t", index=False)
species_out = args.output.replace(".tsv", "_species_level.tsv")
species_matrix.to_csv(species_out, sep="\t", index=False)

print(f"\n===== Feature matrix complete: {DATASET} =====")
print(f"  Timepoint-level rows: {len(feat_matrix)}  → {args.output}")
print(f"  Species-level rows:   {len(species_matrix)}  → {species_out}")
print(f"  Species:              {species_matrix['species'].nunique()}")
print(f"  Features (timepoint): {len(feat_matrix.columns)}")

if not species_matrix.empty:
    print(f"  log2FC range: "
          f"{species_matrix['log2fc'].min():.2f} to "
          f"{species_matrix['log2fc'].max():.2f}")

intarna_status = "included" if args.intarna and os.path.exists(args.intarna or "") \
    else "not included (run IntaRNA and rerun with --intarna to add)"

summary = {
    "dataset":           DATASET,
    "n_timepoint_pairs": len(tp_pairs),
    "n_species":         int(species_matrix["species"].nunique()) if not species_matrix.empty else 0,
    "n_srnas_total":     int(len(snapt_all)),
    "timepoints":        timepoints,
    "pairs":             [f"{a}_vs_{b}" for a, b in tp_pairs],
    "intarna_features":  intarna_status,
    "primary_expression_feature": "FPKM log2FC from SnapT GFF (StringTie)",
}
json_out = args.output.replace(".tsv", "_summary.json")
with open(json_out, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  Summary JSON: {json_out}")

