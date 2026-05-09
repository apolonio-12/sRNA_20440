#!/usr/bin/env python3
"""
sRNA metatranscriptomics pipeline
==================================================
Produces five output tables:
  srna_master.tsv          one row = one sRNA × one sample
  srna_target_pathway.tsv  one row = sRNA → target gene → pathway
  pathway_master.tsv       one row = one pathway × one sample
  srna_presence.tsv        one row = sRNA × sample (presence/absence)
  mechanism_table.tsv      one row = sRNA–pathway pair with consistency flag

Requires: pandas, numpy, scipy, statsmodels
"""

import bz2
import gzip
import os
import re
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ── INPUT PATHS ───────────────────────────────────────────────────────────────

SNAPT_GFF = {
    "2016": "AD/snapt_2016/small_ncRNAs.gff",
    "2017": "AD/snapt_2017/small_ncRNAs.gff",
}
TPM_FILES = {
    "2016": "AD/per_replicate_tpm_2016/tpm_matrix.tsv",
    "2017": "AD/per_replicate_tpm_2017/tpm_matrix.tsv",
}
INTARNA = {
    "2016": "AD/intarna_2016/interactions_significant.tsv",
    "2017": "AD/intarna_2017/interactions_significant.tsv",
}

PRODIGAL_GFF    = "AD/annotation/prodigal.gff"
GENE_TO_UNIREF  = "gene_to_uniref90.tsv"

HUMANN_RNA_GF   = "AD/Total_RNA_genefamilies_humann_table.tsv"
HUMANN_DNA_GF   = "AD/Total_DNA_genefamilies_humann_table.tsv"
HUMANN_RNA_PATH = "AD/Total_RNA_pathabundance_humann_table.tsv"
HUMANN_DNA_PATH = "AD/Total_DNA_pathabundance_humann_table.tsv"

METACYC_RXN_UNIREF = "/home/apolonio/humann_db/utility_mapping/metacyc_reactions_level4ec_only.uniref.bz2"
MAP_PW_NAME        = "/home/apolonio/humann_db/utility_mapping/map_metacyc-pwy_name.txt.gz"
EC_TO_PATHWAY_FILE = "/home/apolonio/humann_db/utility_mapping/metacyc_pathways_structured_filtered_v24_subreactions"

OUT = "FINAL_TABLES"
os.makedirs(OUT, exist_ok=True)

# ── CONDITION ASSIGNMENT ──────────────────────────────────────────────────────

def infer_condition(name: str) -> str:
    n = name.lower()
    if "dry"  in n: return "dry"
    if "wet"  in n: return "wet"
    if "2016" in n: return "dry"
    if "2017" in n: return "wet"
    return "unknown"

def infer_timepoint(name: str) -> str:
    if "2016" in name: return "2016"
    if "2017" in name: return "2017"
    return "unknown"

# ── DIFFERENTIAL EXPRESSION HELPER ───────────────────────────────────────────

def compute_de(
    df: pd.DataFrame,
    fc_threshold: float = 0.5,
    id_col_name: str = "id",
) -> pd.DataFrame:
    """
    Welch t-test + BH-FDR for every row in df (features x samples).
    Samples are split into dry / wet using infer_condition().
    Returns a DataFrame with columns:
        {id_col_name}, log2fc, p_value, padj, direction
    """
    samples  = list(df.columns)
    dry_cols = [c for c in samples if infer_condition(c) == "dry"]
    wet_cols = [c for c in samples if infer_condition(c) == "wet"]

    if not dry_cols or not wet_cols:
        print(
            f"  [WARN] compute_de: could not split samples into dry/wet.\n"
            f"         dry={dry_cols}\n"
            f"         wet={wet_cols}\n"
            f"         Returning NaN for all DE stats."
        )
        return pd.DataFrame({
            id_col_name: df.index,
            "log2fc":    np.nan,
            "p_value":   np.nan,
            "padj":      np.nan,
            "direction": "ns",
        })

    dry_mat = df[dry_cols].values.astype(float)
    wet_mat = df[wet_cols].values.astype(float)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dry_mean = np.nanmean(dry_mat, axis=1)
        wet_mean = np.nanmean(wet_mat, axis=1)

    log2fc = np.log2((wet_mean + 1) / (dry_mean + 1))

    pvals = np.full(len(df), np.nan)
    if len(dry_cols) >= 2 and len(wet_cols) >= 2:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i in range(len(df)):
                _, p = stats.ttest_ind(
                    dry_mat[i], wet_mat[i], equal_var=False, nan_policy="omit"
                )
                pvals[i] = float(p) if not np.ma.is_masked(p) else np.nan

    padj = np.full(len(df), np.nan)
    mask = ~np.isnan(pvals)
    if mask.sum() > 0:
        _, padj[mask], _, _ = multipletests(pvals[mask], method="fdr_bh")

    direction = np.where(
        log2fc >  fc_threshold, "up",
        np.where(log2fc < -fc_threshold, "down", "ns"),
    )

    return pd.DataFrame({
        id_col_name: df.index,
        "log2fc":    log2fc,
        "p_value":   pvals,
        "padj":      padj,
        "direction": direction,
    })

# ── COORDINATE PARSERS ────────────────────────────────────────────────────────

_RE_SRNA_COORDS = re.compile(r"::([^:]+):(\d+)-(\d+)\(([+-])\)")
_RE_MRNA_COORDS = re.compile(r"^([^:]+):(\d+)-(\d+)")

def parse_srna_coords(srna_id: str):
    m = _RE_SRNA_COORDS.search(srna_id)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    return None, None, None, None

def parse_mrna_coords(mrna_id: str):
    m = _RE_MRNA_COORDS.search(mrna_id)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return None, None, None

def coord_key(contig, start, end, strand) -> str:
    """Stable cross-year identifier: contig:start-end:strand"""
    return f"{contig}:{start}-{end}:{strand}" 
    
def overlap_len(a_start, a_end, b_start, b_end):
    return min(a_end, b_end) - max(a_start, b_start)

def same_locus(a, b, min_overlap=10):
    if a["contig"] != b["contig"]:
        return False
    if a["strand"] != b["strand"]:
        return False
    return overlap_len(a["start"], a["end"], b["start"], b["end"]) >= min_overlap

# ── STEP 1 — SNAPT GFF -> sRNA ANNOTATIONS ───────────────────────────────────
# FIX: StringTie re-numbers STRG IDs from 1 each run, so IDs are not
# comparable across years. We build a canonical coord_key and use it as the
# primary identifier throughout. We keep the original srna_id from 2016 where
# available (for readability), falling back to 2017.

_TYPE_MAP = {
    "antisense transcript":  "asRNA",
    "intergenic transcript": "itsRNA",
}

def parse_snapt_gff(path: str, year: str) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            contig, feat_type = p[0], p[2]
            start, end, strand = int(p[3]), int(p[4]), p[6]
            m = re.search(r'gene_id "([^"]+)"', p[8])
            if not m:
                continue
            rows.append({
                "srna_id":  m.group(1),
                "coord_id": coord_key(contig, start, end, strand),
                "type":     _TYPE_MAP.get(feat_type, feat_type),
                "contig":   contig,
                "start":    start,
                "end":      end,
                "strand":   strand,
                "length":   end - start,
                "dataset":  year,
            })
    return pd.DataFrame(rows)

print("Parsing SNAPT GFFs...")
gff_2016 = parse_snapt_gff(SNAPT_GFF["2016"], "2016")
gff_2017 = parse_snapt_gff(SNAPT_GFF["2017"], "2017")

# Build coord_id -> canonical srna_id (prefer 2016, fall back to 2017)
# Combine both years
all_srna = pd.concat([gff_2016, gff_2017], ignore_index=True)

clusters = []   # each cluster = list of rows
cluster_ids = []  # canonical srna_id per cluster

for _, r in all_srna.iterrows():
    assigned = False

    for i, cluster in enumerate(clusters):
        # compare to first member of cluster
        rep = cluster[0]
        if same_locus(r, rep, min_overlap=10):
            cluster.append(r)
            assigned = True
            break

    if not assigned:
        clusters.append([r])
        # prefer 2016 ID if present later
        cluster_ids.append(r["srna_id"])

# Build coord_id -> canonical srna_id
coord_to_srna = {}

for cluster, cid in zip(clusters, cluster_ids):
    # if any member is from 2016, prefer that ID
    ids_2016 = [x["srna_id"] for x in cluster if x["dataset"] == "2016"]
    canonical = ids_2016[0] if ids_2016 else cid

    for x in cluster:
        coord = f'{x["contig"]}:{x["start"]}-{x["end"]}:{x["strand"]}'
        coord_to_srna[coord] = canonical

# Per-year STRG -> coord_id lookups (used to remap TPM rows)
strg2coord = {
    "2016": dict(zip(gff_2016["srna_id"], gff_2016["coord_id"])),
    "2017": dict(zip(gff_2017["srna_id"], gff_2017["coord_id"])),
}

# Canonical annotation table (one row per unique coord_id)
srna_annot = (
    pd.concat([gff_2016, gff_2017], ignore_index=True)
    .sort_values("dataset")           # 2016 rows first so they win drop_duplicates
    .drop_duplicates("coord_id")
    .copy()
)
srna_annot["srna_id"] = srna_annot["coord_id"].map(coord_to_srna)

print(f"  {len(srna_annot)} unique sRNAs (by coordinate)")
print(f"  2016 only: {len(set(gff_2016['coord_id']) - set(gff_2017['coord_id']))}")
print(f"  2017 only: {len(set(gff_2017['coord_id']) - set(gff_2016['coord_id']))}")
print(f"  shared:    {len(set(gff_2016['coord_id']) & set(gff_2017['coord_id']))}")

# Contig-indexed lookup for IntaRNA coordinate matching
srna_by_contig: dict = defaultdict(list)
for _, s in srna_annot.iterrows():
    srna_by_contig[s["contig"]].append(s)

def match_srna_by_coords(contig, start, end, strand):
    best_id, best_overlap = None, 0
    for s in srna_by_contig.get(contig, []):
        overlap = min(s["end"], end) - max(s["start"], start)
        if overlap > best_overlap:
            best_overlap = overlap
            best_id      = s["srna_id"]
    return best_id if best_overlap > 10 else None

# ── STEP 2 — TPM MATRICES -> LONG FORMAT ─────────────────────────────────────
# FIX: remap STRG IDs -> canonical srna_id via coordinates before combining.

print("Loading TPM matrices...")
tpm_long_frames = []
for year, path in TPM_FILES.items():
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "strg_id"})
    df["srna_id"] = df["strg_id"].map(strg2coord[year]).map(coord_to_srna)
    unmapped = df["srna_id"].isna().sum()
    if unmapped:
        print(f"  [WARN] {year}: {unmapped} STRG IDs not found in GFF (dropped)")
    df = df.dropna(subset=["srna_id"])
    long = df.drop(columns=["strg_id"]).melt(
        id_vars="srna_id", var_name="sample_id", value_name="expression"
    )
    long["timepoint"] = year
    long["dataset"]   = year
    long["condition"] = long["sample_id"].apply(infer_condition)
    tpm_long_frames.append(long)

tpm_long = pd.concat(tpm_long_frames, ignore_index=True)
print(f"  {tpm_long['sample_id'].nunique()} samples, "
      f"{tpm_long['srna_id'].nunique()} unique sRNAs (canonical IDs)")

# ── STEP 3 — DIFFERENTIAL EXPRESSION (sRNAs) ─────────────────────────────────
# FIX: pivot uses canonical srna_id so both years' columns are present per row.

print("Computing sRNA differential expression...")
tpm_wide = tpm_long.pivot_table(
    index="srna_id", columns="sample_id", values="expression", aggfunc="mean"
)
print(f"  DE matrix: {tpm_wide.shape[0]} sRNAs x {tpm_wide.shape[1]} samples")
print(f"  Dry cols:  {[c for c in tpm_wide.columns if infer_condition(c)=='dry']}")
print(f"  Wet cols:  {[c for c in tpm_wide.columns if infer_condition(c)=='wet']}")

srna_de = compute_de(tpm_wide, id_col_name="srna_id")
srna_de = srna_de.rename(columns={"direction": "srna_direction"})

n_up   = (srna_de["srna_direction"] == "up").sum()
n_down = (srna_de["srna_direction"] == "down").sum()
print(f"  sRNA DE: {n_up} up, {n_down} down, "
      f"{len(srna_de)-n_up-n_down} ns  (|log2fc|>0.5, no padj filter applied here)")

# ── TABLE 1 — srna_master.tsv ─────────────────────────────────────────────────

print("Building srna_master.tsv...")
srna_master = (
    tpm_long
    .merge(
        srna_annot[["srna_id", "type", "contig", "start", "end", "strand", "length"]],
        on="srna_id", how="left",
    )
    .merge(
        srna_de[["srna_id", "log2fc", "p_value", "padj"]],
        on="srna_id", how="left",
    )
)
srna_master.insert(1, "sequence", np.nan)

srna_master = srna_master[[
    "srna_id", "sequence", "type", "length",
    "sample_id", "dataset", "condition", "timepoint",
    "expression", "log2fc", "p_value", "padj",
    "contig", "start", "end", "strand",
]]
srna_master.to_csv(f"{OUT}/srna_master.tsv", sep="\t", index=False)
print(f"  -> {len(srna_master):,} rows written")

# ── TABLE 2 — srna_presence.tsv ───────────────────────────────────────────────

print("Building srna_presence.tsv...")
presence = tpm_long[["srna_id", "sample_id", "condition"]].copy()
presence["present"] = (tpm_long["expression"] > 0).astype(int)
presence.to_csv(f"{OUT}/srna_presence.tsv", sep="\t", index=False)
print(f"  -> {len(presence):,} rows written")

# ── STEP 4 — PRODIGAL GFF -> GENE INTERVAL INDEX ─────────────────────────────

print("Parsing Prodigal GFF...")
gene_rows = []
with open(PRODIGAL_GFF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        gid  = re.search(r"ID=([^;]+)", p[8])
        name = re.search(r"product=([^;]+)", p[8])
        if not gid:
            continue
        contig       = p[0]
        raw_id       = gid.group(1)
        gene_idx     = raw_id.split("_")[-1]
        full_gene_id = f"{contig}_{gene_idx}"
        gene_rows.append({
            "contig":    contig,
            "start":     int(p[3]),
            "end":       int(p[4]),
            "strand":    p[6],
            "gene_id":   full_gene_id,
            "gene_name": name.group(1) if name else None,
        })

genes_df = pd.DataFrame(gene_rows)
print(f"  {len(genes_df):,} CDS features")

genes_by_contig: dict = defaultdict(list)
for _, g in genes_df.iterrows():
    genes_by_contig[g["contig"]].append(g)

def find_gene(contig, start, end):
    best, best_overlap = None, 0
    for g in genes_by_contig.get(contig, []):
        ov = min(g["end"], end) - max(g["start"], start)
        if ov > best_overlap:
            best_overlap, best = ov, g
    return best

# ── STEP 5 — GENE -> UNIREF MAPPING ──────────────────────────────────────────

print("Loading gene -> UniRef mapping...")
g2u = pd.read_csv(
    GENE_TO_UNIREF,
    sep=r"\s+",
    header=None,
    names=["gene_id", "uniref_raw", "pct_id", "aln_len", "evalue", "score"],
    engine="python",
)
g2u = g2u.sort_values("score", ascending=False).drop_duplicates("gene_id")
g2u["uniref"] = g2u["uniref_raw"].str.split("|").str[0]
g2u = g2u.dropna(subset=["uniref"])
print(f"  {len(g2u):,} genes loaded")
gene2uniref = dict(zip(g2u["gene_id"], g2u["uniref"]))

# ── STEP 6 — BUILD uniref2ec FROM MetaCyc REACTION FILE ──────────────────────

print("Loading UniRef/UniClust -> EC mapping from MetaCyc reactions file...")
uniref2ec = defaultdict(set)
rxn2ec    = {}

with bz2.open(os.path.expanduser(METACYC_RXN_UNIREF), "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        rxn = parts[0]
        ec  = parts[1]
        if not re.match(r"\d+\.\d+\.\d+\.\d+", ec):
            continue
        rxn2ec[rxn] = ec
        for protein_id in parts[2:]:
            uniref2ec[protein_id].add(ec)

print(f"  {len(uniref2ec):,} UniRef/UniClust entries mapped to ECs")
print(f"  {len(rxn2ec):,} reactions have an EC number")

# ── STEP 6.5 — BUILD ec2path VIA REACTION IDs ────────────────────────────────

print("Loading EC -> pathway mapping via reaction IDs...")

_RE_RXN = re.compile(
    r"\b([A-Z0-9][A-Z0-9_-]*RXN[A-Z0-9_-]*"
    r"|RXN[A-Z0-9_-]+)\b"
)
_RE_EC = re.compile(r"\b(\d+\.\d+\.\d+\.\d+)\b")

rxn2path = defaultdict(set)
ec2path  = defaultdict(set)

with open(os.path.expanduser(EC_TO_PATHWAY_FILE)) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        pathway = parts[0].strip()
        rest    = "\t".join(parts[1:])

        for rxn in _RE_RXN.findall(rest):
            rxn2path[rxn].add(pathway)

        for ec in _RE_EC.findall(rest):
            ec2path[ec].add(pathway)

print(f"  {len(rxn2path):,} unique reaction IDs mapped to pathways")

for rxn, ec in rxn2ec.items():
    for pathway in rxn2path.get(rxn, set()):
        ec2path[ec].add(pathway)

print(f"  {len(ec2path):,} EC numbers mapped to pathways (after reaction-ID join)")

# ── STEP 7 — HUMANN GENE FAMILY ABUNDANCES ───────────────────────────────────

print("Loading HUMAnN gene families...")

def load_humann_gf(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df[~df.index.str.contains(r"\|", na=False)]
    return df

rna_gf = load_humann_gf(HUMANN_RNA_GF)
dna_gf = load_humann_gf(HUMANN_DNA_GF)

def gene_mean_abundance(uniref, gf_df: pd.DataFrame) -> float:
    if pd.isna(uniref) or uniref not in gf_df.index:
        return np.nan
    return float(gf_df.loc[uniref].mean())

# ── STEP 8 — HUMANN PATHWAY ABUNDANCES ───────────────────────────────────────
# FIX: HUMAnN pathway IDs are stored as "PWY-123: long name" -- strip suffix.
# FIX: RNA and DNA files use different sample naming (T2 vs 02) so we cannot
#      join per-sample. DNA abundance is stored as the mean across all DNA
#      samples (scalar per pathway) and annotated onto the RNA-based long table.

print("Loading HUMAnN pathways...")

def load_humann_path(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df[~df.index.str.contains(r"\|",                   na=False)]
    df = df[~df.index.str.contains("UNMAPPED|UNINTEGRATED", na=False)]
    df.index = df.index.str.split(":").str[0].str.strip()
    return df

rna_path = load_humann_path(HUMANN_RNA_PATH)
dna_path = load_humann_path(HUMANN_DNA_PATH)

# Scalar mean DNA abundance per pathway
dna_path_mean = dna_path.reindex(rna_path.index).mean(axis=1).rename("abundance_dna")

# Pathway names
pw_names: dict = {}
with gzip.open(MAP_PW_NAME, "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            pw_names[parts[0]] = parts[1]

# Pathway DE
print("Computing pathway differential expression...")
path_de = compute_de(rna_path, id_col_name="pathway_id")
path_de = path_de.rename(columns={"direction": "pathway_direction"})
path_de["pathway_name"]          = path_de["pathway_id"].map(pw_names)
path_de["pathway_abundance_rna"] = rna_path.mean(axis=1).values
path_de["pathway_abundance_dna"] = dna_path_mean.reindex(
    path_de["pathway_id"]
).values

path_de_ids = set(path_de["pathway_id"].values)
path_de_idx = path_de.set_index("pathway_id")

ec2path_pws = set(pw for pws in ec2path.values() for pw in pws)
overlap = ec2path_pws & path_de_ids
print(f"  {len(path_de_ids):,} pathways in HUMAnN table")
print(f"  {len(ec2path_pws):,} pathways reachable via ec2path")
print(f"  {len(overlap):,} pathways overlap (ceiling for Table 3)")

# ── STEP 9 — PARSE INTARNA ───────────────────────────────────────────────────

print("Parsing IntaRNA interactions...")
inter_frames = []
for year, path in INTARNA.items():
    df = pd.read_csv(path, sep="\t")
    df["year"] = year
    inter_frames.append(df)
inter = pd.concat(inter_frames, ignore_index=True)

inter = inter.sort_values("energy").drop_duplicates(subset=["srna_id", "mrna_id"])

pval_col = next(
    (c for c in inter.columns if re.match(r"p.?value", c, re.I)), None
)
inter["interaction_pvalue"] = inter[pval_col] if pval_col else np.nan

coords = inter["srna_id"].apply(
    lambda x: pd.Series(parse_srna_coords(x),
                        index=["srna_contig", "srna_start", "srna_end", "srna_strand"])
)
inter = pd.concat([inter, coords], axis=1)

inter["snapt_srna_id"] = inter.apply(
    lambda r: match_srna_by_coords(
        r["srna_contig"], r["srna_start"], r["srna_end"], r["srna_strand"]
    ),
    axis=1,
)
n_before = len(inter)
inter    = inter.dropna(subset=["snapt_srna_id"])
print(f"  {len(inter):,} / {n_before:,} interactions mapped to a SNAPT sRNA")
if inter.empty:
    raise RuntimeError(
        "No IntaRNA interactions mapped to sRNAs -- check coordinate matching."
    )

target_coords = inter["mrna_id"].apply(
    lambda x: pd.Series(parse_mrna_coords(x),
                        index=["target_contig", "target_start", "target_end"])
)
inter = pd.concat([inter, target_coords], axis=1)

def _map_gene(row):
    g = find_gene(row["target_contig"], row["target_start"], row["target_end"])
    if g is not None:
        return g["gene_id"], g["gene_name"], g["contig"]
    return None, None, None

gene_cols = inter.apply(
    lambda r: pd.Series(
        _map_gene(r),
        index=["target_gene_id", "target_gene_name", "target_contig_gene"]
    ),
    axis=1,
)
inter = pd.concat([inter, gene_cols], axis=1)
inter = inter.dropna(subset=["target_gene_id"])
print(f"  {len(inter):,} interactions after gene mapping")

inter["uniref"] = inter["target_gene_id"].map(gene2uniref)
n_with_uniref   = inter["uniref"].notna().sum()
print(f"  {n_with_uniref:,} / {len(inter):,} interactions have a UniRef ID")

inter["gene_abundance_rna"] = inter["uniref"].apply(
    lambda u: gene_mean_abundance(u, rna_gf)
)
inter["gene_abundance_dna"] = inter["uniref"].apply(
    lambda u: gene_mean_abundance(u, dna_gf)
)

# ── TABLE 3 — srna_target_pathway.tsv ────────────────────────────────────────
# Interactions with no pathway link are retained with NaN pathway fields so
# that the full set of sRNA-target pairs is visible in one place.

print("Building srna_target_pathway.tsv...")

TP_COLS = [
    "srna_id", "target_gene_id", "target_gene_name", "target_contig",
    "interaction_energy", "interaction_pvalue",
    "pathway_id", "pathway_name",
    "gene_abundance_rna", "gene_abundance_dna",
    "pathway_abundance_rna", "pathway_abundance_dna",
    "pathway_log2fc", "pathway_padj", "pathway_direction",
]

_PATHWAY_NAN = {c: np.nan for c in TP_COLS if c not in [
    "srna_id", "target_gene_id", "target_gene_name", "target_contig",
    "interaction_energy", "interaction_pvalue",
    "gene_abundance_rna", "gene_abundance_dna",
]}

tp_rows = []

for _, ir in inter.iterrows():
    base = {
        "srna_id":            ir["snapt_srna_id"],
        "target_gene_id":     ir["target_gene_id"],
        "target_gene_name":   ir["target_gene_name"],
        "target_contig":      ir["target_contig_gene"],
        "interaction_energy": ir["energy"],
        "interaction_pvalue": ir["interaction_pvalue"],
        "gene_abundance_rna": ir["gene_abundance_rna"],
        "gene_abundance_dna": ir["gene_abundance_dna"],
    }

    u = ir["uniref"]
    if pd.isna(u):
        tp_rows.append({**base, **_PATHWAY_NAN})
        continue

    ecs = uniref2ec.get(u, set())
    pathways = set()
    for ec in ecs:
        pathways.update(ec2path.get(ec, set()))
    pathways = pathways & path_de_ids

    if not pathways:
        tp_rows.append({**base, **_PATHWAY_NAN})
        continue

    for pw in pathways:
        pr = path_de_idx.loc[pw]
        tp_rows.append({
            **base,
            "pathway_id":            pw,
            "pathway_name":          pr["pathway_name"],
            "pathway_abundance_rna": pr["pathway_abundance_rna"],
            "pathway_abundance_dna": pr["pathway_abundance_dna"],
            "pathway_log2fc":        pr["log2fc"],
            "pathway_padj":          pr["padj"],
            "pathway_direction":     pr["pathway_direction"],
        })

srna_target_path = pd.DataFrame(tp_rows, columns=TP_COLS)

n_with_pw    = srna_target_path["pathway_id"].notna().sum()
n_without_pw = srna_target_path["pathway_id"].isna().sum()
print(f"  -> {len(srna_target_path):,} rows total "
      f"({n_with_pw:,} with pathway, {n_without_pw:,} without)")

srna_target_path.to_csv(f"{OUT}/srna_target_pathway.tsv", sep="\t", index=False)

# ── TABLE 4 — pathway_master.tsv ─────────────────────────────────────────────

print("Building pathway_master.tsv...")

rna_long = rna_path.reset_index().melt(
    id_vars=rna_path.index.name or "pathway_id",
    var_name="sample_id", value_name="abundance_rna",
)
rna_long.columns = ["pathway_id", "sample_id", "abundance_rna"]
rna_long["condition"]    = rna_long["sample_id"].apply(infer_condition)
rna_long["pathway_name"] = rna_long["pathway_id"].map(pw_names)

# DNA abundance joined per-pathway only (sample names differ between RNA/DNA)
rna_long["abundance_dna"] = rna_long["pathway_id"].map(dna_path_mean)

path_master = rna_long.merge(
    path_de[["pathway_id", "log2fc", "p_value", "padj", "pathway_direction"]],
    on="pathway_id", how="left",
).rename(columns={"pathway_direction": "direction"})

path_master = path_master[[
    "pathway_id", "pathway_name", "sample_id", "condition",
    "abundance_rna", "abundance_dna",
    "log2fc", "p_value", "padj", "direction",
]]
path_master.to_csv(f"{OUT}/pathway_master.tsv", sep="\t", index=False)
print(f"  -> {len(path_master):,} rows written")

# ── TABLE 5 — mechanism_table.tsv ────────────────────────────────────────────

print("Building mechanism_table.tsv...")

MECH_COLS = [
    "srna_id", "srna_log2fc", "srna_direction",
    "pathway_id", "pathway_log2fc", "pathway_direction",
    "consistent",
]

sig = srna_target_path.dropna(subset=["pathway_id"])

if not sig.empty:
    mech = (
        sig[["srna_id", "pathway_id", "pathway_log2fc", "pathway_direction"]]
        .drop_duplicates()
        .merge(
            srna_de[["srna_id", "log2fc", "srna_direction"]].rename(
                columns={"log2fc": "srna_log2fc"}
            ),
            on="srna_id", how="left",
        )
    )
    mech["consistent"] = (
        ((mech["srna_direction"] == "up")   & (mech["pathway_direction"] == "down")) |
        ((mech["srna_direction"] == "down") & (mech["pathway_direction"] == "up"))
    )
    mech = mech[MECH_COLS]
else:
    mech = pd.DataFrame(columns=MECH_COLS)

mech.to_csv(f"{OUT}/mechanism_table.tsv", sep="\t", index=False)
print(f"  -> {len(mech):,} rows written")

# ── SUMMARY ───────────────────────────────────────────────────────────────────

print("\nAll tables written to", OUT)
print(f"   srna_master.tsv         {len(srna_master):>10,} rows")
print(f"   srna_presence.tsv       {len(presence):>10,} rows")
print(f"   srna_target_pathway.tsv {len(srna_target_path):>10,} rows")
print(f"   pathway_master.tsv      {len(path_master):>10,} rows")
print(f"   mechanism_table.tsv     {len(mech):>10,} rows")
