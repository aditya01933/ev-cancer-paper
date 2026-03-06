#!/usr/bin/env python3
"""
Ev Mechanism Analysis: Methylation + TE Experiments
=====================================================
Three experiments to confirm CpG methylation as the Ev mechanism:

  EXP 1: Bisulfite-seq overlay    — Zone3 should have higher CG methylation
  EXP 2: TE annotation overlap    — Zone3 should be enriched in TE windows
  EXP 3: Cross-species CpG O/E   — CpG suppression should predict Zone3%

QUALITY GATES throughout. Each gate is PASS / WARN / FAIL.
CRITICAL gates will abort the script if failed — results would be garbage.

Requirements:
  pip install numpy scipy --break-system-packages

════════════════════════════════════════════════════════════════
INPUT FILES — edit paths below, then run:
  python3 analyze_methylation_te.py
════════════════════════════════════════════════════════════════

  SCAN_JSON    ath_scan_results.json          ← you have this
  ATH_FASTA    ath_chr1.fa                    ← you have this
  MASS_JSON    mass_scan_gcfiltered.json       ← you have this

  METHYL_BED   CG methylation bedGraph        ← download (see below)
  TE_GFF       TAIR10 TE annotation GFF/BED   ← download (see below)
  GENE_GFF     TAIR10 gene annotation GFF     ← download (see below)

  SPECIES_FASTAS  dict of species→fasta path  ← set paths for Exp 3

DOWNLOAD COMMANDS (run in your project folder):
─────────────────────────────────────────────────
  # TE annotation (TAIR10)
  curl "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.gff" -o TAIR10_TE.gff

  # Gene annotation (TAIR10)
  curl "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff" -o TAIR10_genes.gff

  # CG Methylation — Stroud et al. 2013 Col-0 WGBS (GEO: GSE38935)
  # After downloading, extract CG context bedGraph for Chr1
  # Expected format (4 columns, tab-sep):
  #   Chr1  start  end  methylation_fraction
  # where methylation_fraction is 0.0 to 1.0

  # Alternative: use MethylC-seq processed file from TAIR
  wget "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/Methylation/Stroud_2013_CG_methylation_Col0.bedgraph" -O methyl_CG_col0.bedgraph
═══════════════════════════════════════════════════════════════
"""

import json, os, sys, gzip, math
import numpy as np
from scipy.stats import (fisher_exact, pearsonr, spearmanr,
                         mannwhitneyu, kruskal, chi2_contingency)
from collections import defaultdict

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CONFIGURATION — EDIT THESE
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SCAN_JSON  = "ath_scan_results.json"
ATH_FASTA  = "ath_chr1.fa"
MASS_JSON  = "mass_scan_gcfiltered.json"
METHYL_BED = None    # set to bedgraph path if you have one
TE_GFF     = "TAIR10_TE_all.bed"
GENE_GFF   = "pc_genes.gff"   # curl "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff" -o TAIR10_genes.gff

# Species FASTAs — all confirmed present on disk
SPECIES_FASTAS = {
    "Arabidopsis":   "ath_chr1.fa",
    "Rice":          "rice_chr1.fa.gz",
    "Rice_T2T":      "rice_t2t_chr1.fa.gz",
    "Maize":         "maize_chr1.fa.gz",
    "Sorghum":       "sorghum_chr1.fa.gz",
    "Brachypodium":  "brachy_chr1.fa.gz",
    "Soybean":       "soybean_chr1.fa.gz",
    "Physcomitrium": "physco_chr1.fa.gz",
    "Selaginella":   "selaginella.fa.gz",
    "Human":         "human_chr1.fa.gz",
    "Human_T2T":     "human_t2t_chr1.fa.gz",
    "Chlamydomonas": "chlamy_chr1.fa.gz",
    "Drosophila":    "dmel_chr2L.fa.gz",
}

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CONSTANTS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WINDOW_SIZE    = 5000
STEP           = 2500
Z1_T           =  0.7665
Z3_T           = -2.092
GC_MIN         =  0.28
GC_MAX         =  0.60
CHR1_LEN       = 30_427_671
CEN_START      = 14_442_038
CEN_END        = 17_870_129
PERI_START     = 11_000_000
PERI_END       = 21_000_000

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# QUALITY GATE INFRASTRUCTURE
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
_gates = []

def gate(name, condition, value_str, warn_only=False, critical=False):
    """
    Evaluate a quality gate and print result.
    critical=True  → abort script on failure
    warn_only=True → WARN instead of FAIL on failure
    """
    if condition:
        status = "✅ PASS"
    elif warn_only:
        status = "⚠️  WARN"
    else:
        status = "❌ FAIL"
    print(f"    [{status}] {name}")
    print(f"             {value_str}")
    _gates.append({"name": name, "status": status.strip(), "value": value_str})
    if not condition and critical and not warn_only:
        print(f"\n  ══ CRITICAL GATE FAILED — aborting ══")
        print(f"  Fix: {name}\n")
        sys.exit(1)
    return condition

def section(title):
    print(f"\n{'━'*65}")
    print(f"  {title}")
    print(f"{'━'*65}\n")

def subsect(title):
    print(f"\n  ┄ {title} ┄")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# UTILITY FUNCTIONS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def read_fasta_first(path):
    """Return uppercase sequence of first entry in FASTA or gzipped FASTA."""
    opener = gzip.open if path.endswith(".gz") else open
    seq = []
    reading = False
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if reading:
                    break
                reading = True
                print(f"    Header: {line.strip()[:70]}")
                continue
            if reading:
                seq.append(line.strip().upper().replace(" ", ""))
    return "".join(seq)


def gc_content(seq):
    n = len(seq)
    if n == 0:
        return None
    return (seq.count("G") + seq.count("C")) / n


def cpg_oe(seq):
    """CpG observed/expected = (CpG × N) / (C × G)"""
    n = len(seq)
    if n < 100:
        return None
    c = seq.count("C")
    g = seq.count("G")
    if c == 0 or g == 0:
        return None
    cpg = sum(1 for i in range(n - 1) if seq[i:i+2] == "CG")
    return (cpg * n) / (c * g)


def te_fraction_in_window(ws, we, te_ivs_sorted):
    """Fraction of window [ws,we) covered by TE intervals (sorted list)."""
    covered = 0
    for ts, te in te_ivs_sorted:
        if ts >= we:
            break
        ov = min(te, we) - max(ts, ws)
        if ov > 0:
            covered += ov
    return covered / (we - ws)


def gene_fraction_in_window(ws, we, gene_ivs_sorted):
    return te_fraction_in_window(ws, we, gene_ivs_sorted)


def avg_methylation_in_window(ws, we, methyl_sorted):
    """Mean CG methylation fraction in window. methyl_sorted: [(start,end,frac),...]"""
    fracs = [f for s, e, f in methyl_sorted if s < we and e > ws]
    return float(np.mean(fracs)) if fracs else np.nan


def cohens_d(a, b):
    """Cohen's d effect size between two arrays."""
    s = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2)
    return (np.mean(a) - np.mean(b)) / s if s > 0 else 0.0


def woolf_or_ci(a, b, c, d):
    """Odds ratio and 95% CI via Woolf method. Contingency: [[a,b],[c,d]]"""
    if min(a, b, c, d) == 0:
        return None, None, None
    OR = (a * d) / (b * c)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    lo = math.exp(math.log(OR) - 1.96 * se)
    hi = math.exp(math.log(OR) + 1.96 * se)
    return OR, lo, hi


def partial_r(ry_x1, ry_x2, rx1_x2):
    """Partial correlation r(y,x1 | x2)."""
    denom = math.sqrt((1 - ry_x2**2) * (1 - rx1_x2**2))
    if denom == 0:
        return 0.0
    return (ry_x1 - ry_x2 * rx1_x2) / denom


def print_distribution(arr, label, percentiles=(10, 25, 50, 75, 90)):
    valid = arr[~np.isnan(arr)]
    if len(valid) == 0:
        print(f"    {label}: no data")
        return
    pcts = np.percentile(valid, percentiles)
    print(f"    {label}: n={len(valid):6}  mean={np.mean(valid):.4f}"
          f"  sd={np.std(valid):.4f}")
    pct_str = "  ".join(f"P{p}={v:.4f}" for p, v in zip(percentiles, pcts))
    print(f"      Distribution: {pct_str}")


def print_zone_stats(arr, z1m, z2m, z3m, valid_mask, label):
    print(f"\n    {label} by zone:")
    for name, zm in [("Zone1", z1m), ("Zone2", z2m), ("Zone3", z3m)]:
        vals = arr[zm & valid_mask & ~np.isnan(arr)]
        if len(vals) == 0:
            continue
        print(f"      {name}: n={len(vals):5}  mean={np.mean(vals):.4f}"
              f"  med={np.median(vals):.4f}  sd={np.std(vals):.4f}"
              f"  [P25={np.percentile(vals,25):.3f}  P75={np.percentile(vals,75):.3f}]")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# LOAD LOADERS FOR ANNOTATION FILES
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def load_methylation(path):
    """
    Load CG methylation bedGraph for Chr1.
    Accepts: chr start end fraction  (4-col bedGraph)
    Returns sorted list of (start, end, fraction).
    """
    data = []
    opener = gzip.open if path.endswith(".gz") else open
    accepted_chroms = {"1", "Chr1", "chr1", "CHR1", "AT1", "Chromosome1"}
    with opener(path, "rt") as f:
        for i, line in enumerate(f):
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            if parts[0] not in accepted_chroms:
                continue
            try:
                s, e, frac = int(parts[1]), int(parts[2]), float(parts[3])
                if 0.0 <= frac <= 1.0:
                    data.append((s, e, frac))
            except ValueError:
                continue
    data.sort(key=lambda x: x[0])
    return data


def load_te_gff_bed(path):
    """
    Load TE annotation from GFF/GFF3 or BED for Chr1.
    Returns:
      intervals: sorted list of (start, end)
      class_intervals: dict{class_name: sorted[(start,end)]}
    """
    intervals = []
    class_intervals = defaultdict(list)
    accepted_chroms = {"1", "Chr1", "chr1", "CHR1"}
    is_gff = any(path.endswith(x) for x in (".gff", ".gff3", ".gff.gz", ".gff3.gz"))
    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if is_gff:
                if len(parts) < 9:
                    continue
                chrom = parts[0]
                if chrom not in accepted_chroms:
                    continue
                feat = parts[2].lower()
                if not any(k in feat for k in ("transpos", "te", "repeat", "retro", "ltr", "sine", "line", "dna")):
                    continue
                try:
                    s, e = int(parts[3]) - 1, int(parts[4])
                except ValueError:
                    continue
                # Parse class from attributes
                te_class = "Unknown"
                for kv in parts[8].replace(";", " ; ").split(";"):
                    kv = kv.strip()
                    for key in ("Class", "class", "Superfamily", "superfamily", "Family", "family", "Name"):
                        if kv.startswith(key + "=") or kv.startswith(key + " "):
                            te_class = kv.split("=")[-1].strip().split()[0]
                            break
                intervals.append((s, e))
                class_intervals[te_class].append((s, e))
            else:  # BED
                if len(parts) < 3:
                    continue
                if parts[0] not in accepted_chroms:
                    continue
                try:
                    s, e = int(parts[1]), int(parts[2])
                except ValueError:
                    continue
                te_class = parts[3] if len(parts) > 3 else "Unknown"
                intervals.append((s, e))
                class_intervals[te_class].append((s, e))

    intervals.sort(key=lambda x: x[0])
    for cls in class_intervals:
        class_intervals[cls].sort(key=lambda x: x[0])
    return intervals, dict(class_intervals)


def load_genes_gff(path):
    """Load gene intervals from GFF for Chr1. Returns sorted (start, end) list."""
    intervals = []
    accepted_chroms = {"1", "Chr1", "chr1", "CHR1"}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            if parts[0] not in accepted_chroms:
                continue
            if parts[2] == "gene":
                try:
                    intervals.append((int(parts[3]) - 1, int(parts[4])))
                except ValueError:
                    continue
    intervals.sort(key=lambda x: x[0])
    return intervals


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 0: LOAD CORE SCAN DATA
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("STEP 0 — LOAD CORE SCAN DATA")

with open(SCAN_JSON) as f:
    scan = json.load(f)

positions = np.array(scan["positions"])
ev_vals   = np.array(scan["ev_vals"])

gate("Scan loaded — enough windows",
     len(positions) >= 5000,
     f"{len(positions)} windows", critical=True)
gate("Ev range is finite and plausible",
     np.isfinite(ev_vals).all() and ev_vals.min() > -20 and ev_vals.max() < 20,
     f"min={ev_vals.min():.3f}  max={ev_vals.max():.3f}  mean={ev_vals.mean():.3f}",
     critical=True)

z1m = ev_vals > Z1_T
z2m = (ev_vals >= Z3_T) & (ev_vals <= Z1_T)
z3m = ev_vals < Z3_T

gate("Zone proportions in expected range",
     15 < z1m.mean()*100 < 60 and 2 < z3m.mean()*100 < 30,
     f"Zone1={z1m.mean()*100:.1f}%  Zone2={z2m.mean()*100:.1f}%  Zone3={z3m.mean()*100:.1f}%")

# Confirm centromere has high Zone3
cen_mask = (positions >= CEN_START) & (positions <= CEN_END)
arm_mask = ~cen_mask
z3_cen = z3m[cen_mask].mean() * 100
z3_arm = z3m[arm_mask].mean() * 100
gate("Zone3 enriched in centromere vs arms (biological sanity check)",
     z3_cen > z3_arm,
     f"Z3_centromere={z3_cen:.1f}%  Z3_arms={z3_arm:.1f}%")

print(f"\n  Summary:")
print(f"    {len(positions)} total windows   "
      f"Zone1={z1m.sum()} ({z1m.mean()*100:.1f}%)   "
      f"Zone2={z2m.sum()} ({z2m.mean()*100:.1f}%)   "
      f"Zone3={z3m.sum()} ({z3m.mean()*100:.1f}%)")
print(f"    Centromere windows: {cen_mask.sum()}  "
      f"(Z3={z3_cen:.1f}% vs arm Z3={z3_arm:.1f}%)")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STEP 1: SEQUENCE FEATURES FROM FASTA
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("STEP 1 — SEQUENCE FEATURES FROM FASTA")

gc_arr    = np.full(len(positions), np.nan)
cpg_arr   = np.full(len(positions), np.nan)
gc_valid  = np.zeros(len(positions), dtype=bool)
chr1_seq  = None

if not os.path.exists(ATH_FASTA):
    print(f"  FASTA not found: {ATH_FASTA}  — skipping sequence features")
else:
    print(f"  Reading {ATH_FASTA} ...")
    chr1_seq = read_fasta_first(ATH_FASTA)

    gate("Chr1 FASTA length plausible (25–35 Mb)",
         25_000_000 < len(chr1_seq) < 35_000_000,
         f"{len(chr1_seq):,} bp", critical=True)

    ok_chars = set("ACGTN")
    sample_chars = set(chr1_seq[:500_000]) - {"N"}
    gate("FASTA contains only ACGTN",
         sample_chars.issubset(ok_chars),
         f"Extra chars: {sample_chars - ok_chars or 'none'}")

    # Spot-check centromere is AT-rich (known biology)
    cen_seq  = chr1_seq[CEN_START:CEN_START+50_000]
    cen_gc_s = gc_content(cen_seq.replace("N",""))
    gate("Centromere region is AT-rich (GC < 0.43)",
         cen_gc_s < 0.43,
         f"Centromere GC={cen_gc_s:.3f}")

    print(f"\n  Computing GC and CpG O/E for all {len(positions)} windows ...")
    for i, pos in enumerate(positions):
        s = int(pos - WINDOW_SIZE)
        e = int(pos)
        if s < 0 or e > len(chr1_seq):
            continue
        seq = chr1_seq[s:e]
        if seq.count("N") > WINDOW_SIZE * 0.1:
            continue
        gc = gc_content(seq)
        if gc is None:
            continue
        gc_arr[i] = gc
        gc_valid[i] = GC_MIN <= gc <= GC_MAX
        cp = cpg_oe(seq)
        if cp is not None:
            cpg_arr[i] = cp

    n_gcvalid = gc_valid.sum()
    n_cpg     = (~np.isnan(cpg_arr)).sum()
    gate("GC-valid windows > 90% of total",
         n_gcvalid / len(positions) > 0.90,
         f"{n_gcvalid}/{len(positions)} ({n_gcvalid/len(positions)*100:.1f}%)")
    gate("CpG O/E computed for > 90% of windows",
         n_cpg / len(positions) > 0.90,
         f"{n_cpg}/{len(positions)} ({n_cpg/len(positions)*100:.1f}%)")
    gate("CpG O/E values in sane range (0 – 5)",
         np.nanmin(cpg_arr) >= 0 and np.nanmax(cpg_arr) < 5,
         f"min={np.nanmin(cpg_arr):.3f}  max={np.nanmax(cpg_arr):.3f}  mean={np.nanmean(cpg_arr):.3f}")

    # Biological sanity: centromere CpG O/E < arm CpG O/E
    valid_mask = gc_valid & ~np.isnan(cpg_arr)
    cen_cpg = np.nanmean(cpg_arr[cen_mask & valid_mask])
    arm_cpg = np.nanmean(cpg_arr[arm_mask & valid_mask])
    gate("Centromere CpG O/E < arm CpG O/E (methylation sanity)",
         cen_cpg < arm_cpg,
         f"centromere={cen_cpg:.4f}  arm={arm_cpg:.4f}")

    # Replicate previous finding
    r_cpg, p_cpg = pearsonr(ev_vals[valid_mask], cpg_arr[valid_mask])
    rho_cpg, _   = spearmanr(ev_vals[valid_mask], cpg_arr[valid_mask])
    gate("CpG O/E ~ Ev_resid Pearson r > +0.25 (replication of prior result)",
         r_cpg > 0.25,
         f"r={r_cpg:.4f}  rho={rho_cpg:.4f}  p={p_cpg:.2e}")

    subsect("CpG O/E distribution by zone")
    print_zone_stats(cpg_arr, z1m, z2m, z3m, valid_mask, "CpG O/E")

    z1_cpg = cpg_arr[z1m & valid_mask]
    z3_cpg = cpg_arr[z3m & valid_mask]
    d_cpg  = cohens_d(z1_cpg, z3_cpg)
    mw_s, mw_p = mannwhitneyu(z1_cpg, z3_cpg, alternative="two-sided")
    print(f"\n    Z1 vs Z3  Mann-Whitney p={mw_p:.2e}  Cohen's d={d_cpg:.3f}")
    gate("CpG O/E Z1>Z3 difference Cohen's d > 0.30",
         d_cpg > 0.30,
         f"d={d_cpg:.3f}  (>0.3 = medium effect)")

    # Full CpG O/E distribution table (10 decile bins)
    subsect("CpG O/E decile bins → Ev_resid and Zone proportions")
    deciles = np.percentile(cpg_arr[valid_mask], np.arange(0, 110, 10))
    print(f"    {'CpG_O/E bin':22} {'N':>5}  {'Ev_mean':>8}  {'Z1%':>6}  {'Z3%':>6}")
    print("    " + "─"*52)
    for lo, hi in zip(deciles[:-1], deciles[1:]):
        bmask = valid_mask & (cpg_arr >= lo) & (cpg_arr < hi)
        ev_b  = ev_vals[bmask]
        if len(ev_b) < 5:
            continue
        z1p = (ev_b > Z1_T).mean() * 100
        z3p = (ev_b < Z3_T).mean() * 100
        print(f"    [{lo:.3f} – {hi:.3f})         {len(ev_b):>5}  {ev_b.mean():>+8.3f}  {z1p:>5.1f}%  {z3p:>5.1f}%")

    # 1Mb positional bins — Zone1%, Zone3%, CpG O/E mean
    subsect("1Mb bins: Zone1%, Zone3%, CpG O/E mean, Ev mean")
    print(f"    {'Bin':10} {'N':>5}  {'Z1%':>6}  {'Z3%':>6}  {'Ev':>7}  {'CpG_O/E':>9}  {'region'}")
    print("    " + "─"*65)
    for b in range(0, 31_000_000, 1_000_000):
        bmask = (positions >= b) & (positions < b + 1_000_000) & valid_mask
        ev_b  = ev_vals[bmask]
        cp_b  = cpg_arr[bmask]
        if len(ev_b) < 5:
            continue
        region = ""
        if b < CEN_END and b + 1_000_000 > CEN_START:
            region = "◀ CENTROMERE"
        elif b < PERI_END and b + 1_000_000 > PERI_START:
            region = "◁ pericen"
        z1p = (ev_b > Z1_T).mean() * 100
        z3p = (ev_b < Z3_T).mean() * 100
        print(f"    {b//1e6:.0f}–{(b+1e6)//1e6:.0f}Mb      {len(ev_b):>5}  "
              f"{z1p:>5.1f}%  {z3p:>5.1f}%  {ev_b.mean():>+7.3f}  "
              f"{np.nanmean(cp_b):>9.4f}  {region}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# EXPERIMENT 1: BISULFITE-SEQ METHYLATION
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("EXPERIMENT 1 — CG Methylation Overlap")

win_methyl = np.full(len(positions), np.nan)

if not METHYL_BED or not os.path.exists(METHYL_BED):
    print(f"""  METHYL_BED not set or not found.

  Download options:
  ─────────────────
  # Option A: TAIR-hosted processed file (if available)
  wget "https://www.arabidopsis.org/.../Stroud_2013_CG_methylation_Col0.bedgraph"

  # Option B: Process from GEO GSE38935 (Stroud et al. 2013)
  #   Download SRR files, trim, align to TAIR10, extract CG context with Bismark
  #   Output bedGraph: Chr1  start  end  CG_fraction

  # Option C: Use any Col-0 WGBS CG bedGraph (e.g. from SRA or ArrayExpress)

  Expected format (tab-separated, Chr1 only):
    Chr1  100  101  0.95
    Chr1  150  151  0.02
    ...

  Set METHYL_BED = "your_file.bedgraph" and rerun.
""")
else:
    print(f"  Loading {METHYL_BED} ...")
    methyl_data = load_methylation(METHYL_BED)

    gate("Methylation file loaded",
         len(methyl_data) >= 1000,
         f"{len(methyl_data):,} CG sites", critical=True)
    gate("CG fraction values in [0,1]",
         all(0.0 <= x[2] <= 1.0 for x in methyl_data[:1000]),
         f"first-1000 fractions OK")
    gate("Methylation data spans Chr1",
         methyl_data[-1][0] > 20_000_000,
         f"last position: {methyl_data[-1][0]:,}")

    fracs = np.array([x[2] for x in methyl_data])
    print(f"\n  Genome-wide CG methylation stats:")
    print_distribution(fracs, "CG fraction")
    pct_high = (fracs > 0.7).mean() * 100
    pct_low  = (fracs < 0.1).mean() * 100
    print(f"    Highly methylated (>0.7): {pct_high:.1f}%")
    print(f"    Unmethylated (<0.1):      {pct_low:.1f}%")
    gate("Bimodal methylation expected: >20% low AND >20% high",
         pct_low > 20 and pct_high > 20,
         f"low={pct_low:.1f}%  high={pct_high:.1f}%")

    # Biological sanity: centromere should be highly methylated
    cen_sites  = [x[2] for x in methyl_data if CEN_START <= x[0] <= CEN_END]
    arm_sites  = [x[2] for x in methyl_data if x[0] < CEN_START or x[0] > CEN_END]
    if cen_sites:
        gate("Centromere CG methylation > arm methylation (sanity)",
             np.mean(cen_sites) > np.mean(arm_sites),
             f"centromere={np.mean(cen_sites):.3f}  arm={np.mean(arm_sites):.3f}")

    # Per-window average methylation
    print(f"\n  Computing per-window mean methylation ({len(positions)} windows) ...")
    for i, pos in enumerate(positions):
        ws, we = int(pos - WINDOW_SIZE), int(pos)
        win_methyl[i] = avg_methylation_in_window(ws, we, methyl_data)
        if i % 3000 == 0:
            print(f"    {i}/{len(positions)} ...")

    methyl_valid = ~np.isnan(win_methyl)
    gate("Methylation coverage > 80% of windows",
         methyl_valid.mean() > 0.80,
         f"{methyl_valid.mean()*100:.1f}% windows have methylation data")

    # Full distribution by zone
    subsect("Per-window methylation by zone")
    print_zone_stats(win_methyl, z1m, z2m, z3m, methyl_valid & gc_valid, "CG methylation")

    z1_m = win_methyl[z1m & methyl_valid & gc_valid]
    z2_m = win_methyl[z2m & methyl_valid & gc_valid]
    z3_m = win_methyl[z3m & methyl_valid & gc_valid]

    d_meth = cohens_d(z3_m, z1_m)
    _, p_mw = mannwhitneyu(z3_m, z1_m, alternative="greater")
    print(f"\n    Zone3 > Zone1 methylation: p={p_mw:.2e}  Cohen's d={d_meth:.3f}")

    gate("Zone3 higher methylation than Zone1 (p < 0.001)",
         p_mw < 0.001,
         f"p={p_mw:.2e}")
    gate("Effect size meaningful (Cohen's d > 0.30)",
         d_meth > 0.30,
         f"d={d_meth:.3f}")

    # Kruskal across all three zones
    kw_stat, kw_p = kruskal(z1_m, z2_m, z3_m)
    gate("Kruskal-Wallis across zones significant (p < 1e-10)",
         kw_p < 1e-10,
         f"H={kw_stat:.2f}  p={kw_p:.2e}")

    # Pearson + Spearman: Ev_resid ~ methylation
    mask = methyl_valid & gc_valid
    r_m,   p_r  = pearsonr(ev_vals[mask], win_methyl[mask])
    rho_m, p_rho= spearmanr(ev_vals[mask], win_methyl[mask])
    print(f"\n    Ev_resid ~ CG methylation:")
    print(f"      Pearson  r={r_m:+.4f}  p={p_r:.2e}")
    print(f"      Spearman r={rho_m:+.4f}  p={p_rho:.2e}")
    gate("Ev_resid negatively correlated with methylation (r < -0.25)",
         r_m < -0.25,
         f"r={r_m:.4f}")

    # Methylation quartile → Ev_resid: should be strictly decreasing
    subsect("Methylation quartile bins → Ev_resid (must be monotonic)")
    q25, q50, q75 = np.percentile(win_methyl[methyl_valid], [25, 50, 75])
    qs = [(f"Q1 low  meth ≤{q25:.2f}", 0.0,  q25),
          (f"Q2      meth ≤{q50:.2f}", q25,  q50),
          (f"Q3      meth ≤{q75:.2f}", q50,  q75),
          (f"Q4 high meth  >{q75:.2f}", q75, 1.01)]
    prev, monotone = None, True
    qrows = []
    for qlabel, qlo, qhi in qs:
        qmask = methyl_valid & gc_valid & (win_methyl >= qlo) & (win_methyl < qhi)
        ev_q  = ev_vals[qmask]
        if len(ev_q) == 0:
            continue
        m_ev  = ev_q.mean()
        z1p   = (ev_q > Z1_T).mean() * 100
        z3p   = (ev_q < Z3_T).mean() * 100
        print(f"    {qlabel:28}: N={len(ev_q):5}  Ev={m_ev:+.3f}"
              f"  Z1={z1p:.1f}%  Z3={z3p:.1f}%")
        qrows.append(m_ev)
        if prev is not None and m_ev > prev:
            monotone = False
        prev = m_ev
    gate("Ev_resid monotonically decreases across methylation quartiles",
         monotone,
         "monotonic ✓" if monotone else "NOT monotonic — unexpected pattern")

    # Methylation decile bins for finer resolution
    subsect("Methylation decile bins → Zone proportions")
    deciles = np.percentile(win_methyl[methyl_valid], np.arange(0, 110, 10))
    print(f"    {'Methylation bin':24} {'N':>5}  {'Ev_mean':>8}  {'Z1%':>6}  {'Z3%':>6}")
    print("    " + "─"*55)
    for lo, hi in zip(deciles[:-1], deciles[1:]):
        bmask = methyl_valid & gc_valid & (win_methyl >= lo) & (win_methyl < hi)
        ev_b  = ev_vals[bmask]
        if len(ev_b) < 5:
            continue
        z1p = (ev_b > Z1_T).mean() * 100
        z3p = (ev_b < Z3_T).mean() * 100
        print(f"    [{lo:.3f} – {hi:.3f})          {len(ev_b):>5}  "
              f"{ev_b.mean():>+8.3f}  {z1p:>5.1f}%  {z3p:>5.1f}%")

    # Combine CpG O/E + methylation: are they correlated?
    if chr1_seq is not None:
        combined_mask = methyl_valid & gc_valid & ~np.isnan(cpg_arr)
        r_cpg_meth, _ = pearsonr(cpg_arr[combined_mask], win_methyl[combined_mask])
        print(f"\n    r(CpG O/E, methylation) = {r_cpg_meth:+.4f}")
        gate("CpG O/E negatively correlated with methylation (r < -0.20)",
             r_cpg_meth < -0.20,
             f"r={r_cpg_meth:.4f}  (expected: high CpG O/E = unmethylated)")

        # Partial correlations: CpG O/E and methylation both predicting Ev
        r_ev_cpg, _ = pearsonr(ev_vals[combined_mask], cpg_arr[combined_mask])
        r_ev_meth,_ = pearsonr(ev_vals[combined_mask], win_methyl[combined_mask])
        pr_cpg  = partial_r(r_ev_cpg, r_ev_meth, r_cpg_meth)
        pr_meth = partial_r(r_ev_meth, r_ev_cpg, r_cpg_meth)
        print(f"\n    Partial correlations with Ev_resid:")
        print(f"      r(Ev, CpG_O/E | methylation) = {pr_cpg:+.4f}")
        print(f"      r(Ev, methylation | CpG_O/E)  = {pr_meth:+.4f}")
        gate("CpG O/E has independent partial r > 0.10 (not fully redundant with methylation)",
             abs(pr_cpg) > 0.10,
             f"partial_r={pr_cpg:.4f}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# EXPERIMENT 2: TE ANNOTATION OVERLAP
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("EXPERIMENT 2 — TE Annotation Overlap")

te_frac_arr  = np.zeros(len(positions))
te_intervals = []
te_classes   = {}

if not TE_GFF or not os.path.exists(TE_GFF):
    print(f"""  TE_GFF not set or not found.

  Download TAIR10 TE annotation:
    curl "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.gff" -o TAIR10_TE.gff
  or
    wget "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_TE.bed" -O TAIR10_TE.bed

  Then set TE_GFF = "TAIR10_TE.gff" and rerun.
""")
else:
    print(f"  Loading {TE_GFF} ...")
    te_intervals, te_classes = load_te_gff_bed(TE_GFF)

    gate("TE annotation loaded",
         len(te_intervals) >= 200,
         f"{len(te_intervals):,} TE intervals on Chr1", critical=True)

    total_te_bp = sum(e - s for s, e in te_intervals)
    te_cov_pct  = total_te_bp / CHR1_LEN * 100
    gate("TE coverage plausible for Arabidopsis Chr1 (5–25%)",
         5 < te_cov_pct < 25,
         f"{te_cov_pct:.1f}% of Chr1 covered")

    n_cen_tes = sum(1 for s, e in te_intervals if s >= CEN_START and e <= CEN_END)
    gate("TEs present in centromere (biological sanity)",
         n_cen_tes >= 5,
         f"{n_cen_tes} TEs in centromere core")

    print(f"\n  TE class breakdown (top 15 by count):")
    cls_counts = sorted(te_classes.items(), key=lambda x: -len(x[1]))
    for cls, ivs in cls_counts[:15]:
        bp = sum(e - s for s, e in ivs)
        print(f"    {cls:30}: n={len(ivs):5}  coverage={bp/1e6:.2f} Mb")

    # Compute per-window TE fraction
    print(f"\n  Computing TE fraction per window ({len(positions)} windows) ...")
    for i, pos in enumerate(positions):
        ws, we = int(pos - WINDOW_SIZE), int(pos)
        te_frac_arr[i] = te_fraction_in_window(ws, we, te_intervals)
        if i % 3000 == 0:
            print(f"    {i}/{len(positions)} ...")

    # Full distribution
    subsect("TE fraction distribution by zone")
    print_zone_stats(te_frac_arr, z1m, z2m, z3m, gc_valid, "TE fraction")

    z1_te = te_frac_arr[z1m & gc_valid]
    z2_te = te_frac_arr[z2m & gc_valid]
    z3_te = te_frac_arr[z3m & gc_valid]
    d_te  = cohens_d(z3_te, z1_te)
    _, p_te_mw = mannwhitneyu(z3_te, z1_te, alternative="greater")
    kw_s, kw_p = kruskal(z1_te, z2_te, z3_te)
    print(f"\n    Zone3 > Zone1 TE fraction: p={p_te_mw:.2e}  Cohen's d={d_te:.3f}")
    print(f"    Kruskal-Wallis 3 zones:    p={kw_p:.2e}")

    gate("Zone3 higher TE fraction than Zone1 (p < 0.001)",
         p_te_mw < 0.001, f"p={p_te_mw:.2e}")
    gate("TE fraction effect size meaningful (d > 0.30)",
         d_te > 0.30, f"Cohen's d={d_te:.3f}")

    # Pearson + Spearman
    r_te,   _ = pearsonr(ev_vals[gc_valid], te_frac_arr[gc_valid])
    rho_te, _ = spearmanr(ev_vals[gc_valid], te_frac_arr[gc_valid])
    print(f"\n    Ev_resid ~ TE fraction:  r={r_te:+.4f}  rho={rho_te:+.4f}")
    gate("Ev_resid negatively correlated with TE fraction (r < -0.20)",
         r_te < -0.20, f"r={r_te:.4f}")

    # TE fraction decile bins
    subsect("TE fraction decile bins → Zone proportions")
    te_dec = np.percentile(te_frac_arr[gc_valid], np.arange(0, 110, 10))
    print(f"    {'TE fraction bin':22} {'N':>5}  {'Ev_mean':>8}  {'Z1%':>6}  {'Z3%':>6}")
    print("    " + "─"*52)
    for lo, hi in zip(te_dec[:-1], te_dec[1:]):
        bmask = gc_valid & (te_frac_arr >= lo) & (te_frac_arr < hi)
        ev_b  = ev_vals[bmask]
        if len(ev_b) < 5:
            continue
        z1p = (ev_b > Z1_T).mean() * 100
        z3p = (ev_b < Z3_T).mean() * 100
        print(f"    [{lo:.2f} – {hi:.2f})            {len(ev_b):>5}  "
              f"{ev_b.mean():>+8.3f}  {z1p:>5.1f}%  {z3p:>5.1f}%")

    # Fisher exact: Zone3 enrichment in TE-high (>50%) vs TE-low (<10%)
    te_hi_mask = gc_valid & (te_frac_arr >= 0.50)
    te_lo_mask = gc_valid & (te_frac_arr < 0.10)
    for thresh_label, hi_mask, lo_mask in [
        ("TE≥50% vs TE<10%", te_hi_mask, te_lo_mask),
        ("TE≥30% vs TE<30%",
         gc_valid & (te_frac_arr >= 0.30),
         gc_valid & (te_frac_arr < 0.30)),
    ]:
        a = (z3m & hi_mask).sum()
        b = (~z3m & hi_mask).sum()
        c = (z3m & lo_mask).sum()
        d = (~z3m & lo_mask).sum()
        if min(a, b, c, d) > 0:
            OR_v, ci_lo, ci_hi = woolf_or_ci(a, b, c, d)
            _, p_fe = fisher_exact([[a, b], [c, d]])
            print(f"\n    Zone3 enrichment {thresh_label}:")
            print(f"      a={a} b={b} c={c} d={d}")
            print(f"      OR={OR_v:.3f}  95%CI [{ci_lo:.2f}–{ci_hi:.2f}]  p={p_fe:.2e}")
            gate(f"Zone3 enriched in TE-dense windows ({thresh_label}, OR>2)",
                 OR_v > 2, f"OR={OR_v:.3f}  p={p_fe:.2e}")

    # TE class-by-class correlation with Ev_resid
    subsect("Per TE class: r(Ev_resid, TE_class_fraction)")
    print(f"    {'TE class':32} {'n_TEs':>6}  {'r_Ev':>7}  {'rho_Ev':>8}  {'Z3_mean':>8}")
    print("    " + "─"*65)
    cls_results = []
    for cls, ivs in cls_counts[:20]:
        cls_frac = np.zeros(len(positions))
        for i, pos in enumerate(positions):
            ws, we = int(pos - WINDOW_SIZE), int(pos)
            cls_frac[i] = te_fraction_in_window(ws, we, ivs)
        r_cls,_   = pearsonr(ev_vals[gc_valid], cls_frac[gc_valid])
        rho_cls,_ = spearmanr(ev_vals[gc_valid], cls_frac[gc_valid])
        z3_cls    = cls_frac[z3m & gc_valid].mean()
        print(f"    {cls:32} {len(ivs):>6}  {r_cls:>+7.4f}  {rho_cls:>+8.4f}  {z3_cls:>8.4f}")
        cls_results.append((cls, r_cls, rho_cls))

    # Which class has strongest negative r with Ev?
    best_cls = min(cls_results, key=lambda x: x[1])
    print(f"\n    Strongest negative r with Ev: {best_cls[0]} (r={best_cls[1]:.4f})")

    # If methylation data available: TE fraction vs methylation
    if np.any(~np.isnan(win_methyl)):
        combined = gc_valid & methyl_valid & (te_frac_arr > 0)
        if combined.sum() > 100:
            r_te_meth, _ = pearsonr(te_frac_arr[combined], win_methyl[combined])
            print(f"\n    r(TE fraction, CG methylation) = {r_te_meth:+.4f}")
            gate("TE fraction positively correlated with methylation (r > 0.20)",
                 r_te_meth > 0.20,
                 f"r={r_te_meth:.4f}  (TEs should be methylated)")

    # Centromere vs arm TE analysis
    subsect("TE fraction: Centromere vs pericentromere vs arms")
    for label, mask_ in [
        ("Centromere core", cen_mask & gc_valid),
        ("Pericentromere",  (positions >= PERI_START) & (positions < PERI_END)
                           & ~cen_mask & gc_valid),
        ("Arms",           ((positions < PERI_START) | (positions >= PERI_END)) & gc_valid),
    ]:
        tf = te_frac_arr[mask_]
        ev = ev_vals[mask_]
        z3p = z3m[mask_].mean() * 100
        print(f"    {label:25}: N={len(tf):5}  TE_frac={tf.mean():.3f}"
              f"  Z3={z3p:.1f}%  Ev_mean={ev.mean():+.3f}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# BONUS: GENE DENSITY
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("BONUS — Gene Density vs Zone1%")

gene_frac_arr = np.zeros(len(positions))

if not GENE_GFF or not os.path.exists(GENE_GFF):
    print(f"""  GENE_GFF not set or not found.

  Download:
    curl "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff" -o TAIR10_genes.gff

  Set GENE_GFF = "TAIR10_genes.gff" and rerun.
""")
else:
    print(f"  Loading {GENE_GFF} ...")
    gene_ivs = load_genes_gff(GENE_GFF)

    gate("Gene annotation loaded",
         len(gene_ivs) >= 500,
         f"{len(gene_ivs):,} genes on Chr1", critical=False)

    gene_cov_pct = sum(e-s for s,e in gene_ivs) / CHR1_LEN * 100
    gate("Gene coverage plausible (25–55% for Arabidopsis)",
         25 < gene_cov_pct < 55,
         f"{gene_cov_pct:.1f}%")

    cen_genes = sum(1 for s,e in gene_ivs if s >= CEN_START and e <= CEN_END)
    gate("Centromere is gene-poor (< 20 genes in core)",
         cen_genes < 20,
         f"{cen_genes} genes in centromere core")

    print(f"  Computing gene fraction per window ...")
    for i, pos in enumerate(positions):
        ws, we = int(pos - WINDOW_SIZE), int(pos)
        gene_frac_arr[i] = gene_fraction_in_window(ws, we, gene_ivs)
        if i % 3000 == 0:
            print(f"    {i}/{len(positions)} ...")

    r_gene,   _ = pearsonr(ev_vals[gc_valid], gene_frac_arr[gc_valid])
    rho_gene, _ = spearmanr(ev_vals[gc_valid], gene_frac_arr[gc_valid])
    print(f"\n  Ev_resid ~ gene fraction: r={r_gene:+.4f}  rho={rho_gene:+.4f}")

    subsect("Gene fraction by zone")
    print_zone_stats(gene_frac_arr, z1m, z2m, z3m, gc_valid, "Gene fraction")

    gate("Zone1 enriched in gene-dense windows (r > 0.10)",
         r_gene > 0.10,
         f"r={r_gene:.4f}")

    # Gene fraction decile bins
    subsect("Gene fraction decile bins → Zone proportions")
    gene_dec = np.percentile(gene_frac_arr[gc_valid], np.arange(0, 110, 10))
    print(f"    {'Gene fraction bin':22} {'N':>5}  {'Ev_mean':>8}  {'Z1%':>6}  {'Z3%':>6}")
    print("    " + "─"*52)
    for lo, hi in zip(gene_dec[:-1], gene_dec[1:]):
        bmask = gc_valid & (gene_frac_arr >= lo) & (gene_frac_arr < hi)
        ev_b  = ev_vals[bmask]
        if len(ev_b) < 5:
            continue
        z1p = (ev_b > Z1_T).mean() * 100
        z3p = (ev_b < Z3_T).mean() * 100
        print(f"    [{lo:.2f} – {hi:.2f})            {len(ev_b):>5}  "
              f"{ev_b.mean():>+8.3f}  {z1p:>5.1f}%  {z3p:>5.1f}%")

    # If TE data also present: variance partitioning
    if len(te_intervals) > 0 and chr1_seq is not None:
        subsect("Variance partitioning: how much of Ev is explained by each predictor?")
        mask_all = gc_valid & ~np.isnan(cpg_arr)
        ev_s  = ev_vals[mask_all]
        cpg_s = cpg_arr[mask_all]
        te_s  = te_frac_arr[mask_all]
        gn_s  = gene_frac_arr[mask_all]
        gc_s  = gc_arr[mask_all]

        print(f"    Predictor        r         r²      variance_explained")
        print("    " + "─"*52)
        for pred_name, pred in [("GC content",    gc_s),
                                  ("CpG O/E",       cpg_s),
                                  ("TE fraction",   te_s),
                                  ("Gene fraction", gn_s)]:
            r_v, _ = pearsonr(ev_s, pred)
            print(f"    {pred_name:16} {r_v:>+8.4f}  {r_v**2:>6.4f}  {r_v**2*100:.1f}%")

        print(f"\n    Note: predictors are correlated — total ≠ sum of r²")
        print(f"    True multivariate R² requires regression (beyond scope here)")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# EXPERIMENT 3: CROSS-SPECIES CpG O/E vs Zone3%
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("EXPERIMENT 3 — Cross-Species CpG O/E vs Zone3%")

# Load mass scan data
sp_ev_stats = {}
if os.path.exists(MASS_JSON):
    with open(MASS_JSON) as f:
        sp_ev_stats = json.load(f)
    print(f"  Mass scan loaded: {len(sp_ev_stats)} species")
    gate("Mass scan has ≥ 7 species",
         len(sp_ev_stats) >= 7,
         f"{len(sp_ev_stats)} species")
else:
    print(f"  Mass scan not found: {MASS_JSON}")

# Compute CpG O/E from species FASTAs
print(f"\n  Computing CpG O/E per species (500-window sample) ...\n")
sp_cpg_results = {}

for sp, fasta_path in SPECIES_FASTAS.items():
    if not fasta_path or not os.path.exists(fasta_path):
        print(f"  {sp:20}: FASTA not provided / not found")
        continue
    try:
        seq = read_fasta_first(fasta_path)
    except Exception as ex:
        print(f"  {sp:20}: ERROR — {ex}")
        continue
    if len(seq) < 500_000:
        print(f"  {sp:20}: sequence too short ({len(seq):,}) — skip")
        continue

    # Sample up to 500 GC-valid windows
    step_s = max(WINDOW_SIZE, len(seq) // 500)
    cpgs, gcs_list = [], []
    n_skip_gc, n_skip_n = 0, 0
    for start in range(0, len(seq) - WINDOW_SIZE, step_s):
        w = seq[start : start + WINDOW_SIZE]
        if w.count("N") > WINDOW_SIZE * 0.1:
            n_skip_n += 1
            continue
        gc = gc_content(w)
        if gc is None or not (GC_MIN <= gc <= GC_MAX):
            n_skip_gc += 1
            continue
        cp = cpg_oe(w)
        if cp is not None:
            cpgs.append(cp)
            gcs_list.append(gc)

    if len(cpgs) < 30:
        print(f"  {sp:20}: too few valid windows ({len(cpgs)}) — skip")
        continue

    sp_cpg_results[sp] = {
        "mean_cpg_oe":   float(np.mean(cpgs)),
        "median_cpg_oe": float(np.median(cpgs)),
        "sd_cpg_oe":     float(np.std(cpgs)),
        "mean_gc":       float(np.mean(gcs_list)),
        "n_windows":     len(cpgs),
        "n_skip_gc":     n_skip_gc,
        "n_skip_n":      n_skip_n,
    }
    print(f"  {sp:20}: CpG_O/E={np.mean(cpgs):.4f}±{np.std(cpgs):.4f}"
          f"  GC={np.mean(gcs_list):.3f}  N={len(cpgs)}"
          f"  (skip_gc={n_skip_gc} skip_N={n_skip_n})")

gate("CpG O/E computed for ≥ 4 species",
     len(sp_cpg_results) >= 4,
     f"{len(sp_cpg_results)} species computed",
     warn_only=True)

# Cross-tabulation table
print(f"\n  Cross-species table:")
print(f"  {'Species':16} {'Clade':12} {'CpG_O/E':>10} {'GC':>6}  {'Zone3%':>8} {'Z3_mean':>9} {'Z1%':>7}")
print("  " + "─"*75)

xs_cpg, xs_z3pct, xs_names, xs_clades = [], [], [], []

for sp, cpg_info in sorted(sp_cpg_results.items()):
    clade   = sp_ev_stats.get(sp, {}).get("clade", "Unknown")
    z3_pct  = sp_ev_stats.get(sp, {}).get("z3_pct",  float("nan"))
    z3_mean = sp_ev_stats.get(sp, {}).get("z3_mean", float("nan"))
    z1_pct  = sp_ev_stats.get(sp, {}).get("z1_pct",  float("nan"))
    cpg_v   = cpg_info["mean_cpg_oe"]
    gc_v    = cpg_info["mean_gc"]
    print(f"  {sp:16} {clade:12} {cpg_v:>10.4f} {gc_v:>6.3f}  {z3_pct:>7.2f}%"
          f" {z3_mean:>9.3f} {z1_pct:>6.2f}%")
    if not math.isnan(z3_pct):
        xs_cpg.append(cpg_v)
        xs_z3pct.append(z3_pct)
        xs_names.append(sp)
        xs_clades.append(clade)

if len(xs_cpg) >= 4:
    xs_cpg   = np.array(xs_cpg)
    xs_z3pct = np.array(xs_z3pct)

    r_xs,   p_xs   = pearsonr(xs_cpg, xs_z3pct)
    rho_xs, p_rho  = spearmanr(xs_cpg, xs_z3pct)
    print(f"\n  Cross-species Pearson  r(CpG_O/E, Zone3%) = {r_xs:+.4f}  p={p_xs:.3f}")
    print(f"  Cross-species Spearman r(CpG_O/E, Zone3%) = {rho_xs:+.4f}  p={p_rho:.3f}")
    print(f"  Expected sign: NEGATIVE (lower CpG O/E = more CpG suppression = more TEs = higher Zone3%)")

    gate("Cross-species: CpG O/E negatively predicts Zone3% (r < -0.30)",
         r_xs < -0.30,
         f"r={r_xs:.4f}  rho={rho_xs:.4f}  p={p_xs:.3f}",
         warn_only=len(xs_cpg) < 6)   # warn only if few species

    # Monocot vs dicot CpG O/E comparison
    mono_cpg  = [xs_cpg[i] for i, c in enumerate(xs_clades) if c == "Monocot"]
    dicot_cpg = [xs_cpg[i] for i, c in enumerate(xs_clades) if c == "Dicot"]
    mono_z3   = [xs_z3pct[i] for i, c in enumerate(xs_clades) if c == "Monocot"]
    dicot_z3  = [xs_z3pct[i] for i, c in enumerate(xs_clades) if c == "Dicot"]

    if mono_cpg and dicot_cpg:
        print(f"\n  Monocots: CpG_O/E={np.mean(mono_cpg):.4f}  Zone3%={np.mean(mono_z3):.2f}%  n={len(mono_cpg)}")
        print(f"  Dicots:   CpG_O/E={np.mean(dicot_cpg):.4f}  Zone3%={np.mean(dicot_z3):.2f}%  n={len(dicot_cpg)}")
        gate("Monocots have lower CpG O/E than dicots (consistent with higher Zone3%)",
             np.mean(mono_cpg) < np.mean(dicot_cpg),
             f"mono={np.mean(mono_cpg):.4f}  dicot={np.mean(dicot_cpg):.4f}",
             warn_only=True)

    # Per-species CpG O/E distribution details
    subsect("Per-species CpG O/E distribution")
    print(f"  {'Species':16} {'mean':>8} {'median':>8} {'sd':>7} {'n':>5}")
    print("  " + "─"*48)
    for sp, info in sorted(sp_cpg_results.items()):
        print(f"  {sp:16} {info['mean_cpg_oe']:>8.4f} {info['median_cpg_oe']:>8.4f}"
              f" {info['sd_cpg_oe']:>7.4f} {info['n_windows']:>5}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# QUALITY GATE SUMMARY
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("QUALITY GATE SUMMARY")

n_pass  = sum(1 for g in _gates if "PASS" in g["status"])
n_warn  = sum(1 for g in _gates if "WARN" in g["status"])
n_fail  = sum(1 for g in _gates if "FAIL" in g["status"])
total   = len(_gates)

print(f"  Total gates evaluated: {total}")
print(f"  ✅ PASS:   {n_pass}")
print(f"  ⚠️  WARN:   {n_warn}")
print(f"  ❌ FAIL:   {n_fail}")

if n_fail:
    print(f"\n  ❌ FAILED GATES:")
    for g in _gates:
        if "FAIL" in g["status"]:
            print(f"     • {g['name']}")
            print(f"       {g['value']}")
if n_warn:
    print(f"\n  ⚠️  WARNINGS:")
    for g in _gates:
        if "WARN" in g["status"]:
            print(f"     • {g['name']}")
            print(f"       {g['value']}")

if n_fail == 0:
    print(f"\n  ✅ No failures — results are trustworthy.")
else:
    print(f"\n  ❌ {n_fail} gate(s) failed — review those sections before drawing conclusions.")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SAVE ALL RESULTS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("SAVING RESULTS → methylation_te_results.json")

def safe(v):
    if isinstance(v, (np.floating, np.integer)):
        return float(v)
    if isinstance(v, np.ndarray):
        return v.tolist()
    return v

output = {
    "quality_gates":  _gates,
    "gate_summary":   {"pass": n_pass, "warn": n_warn, "fail": n_fail, "total": total},
    "zone_counts":    {"Z1": int(z1m.sum()), "Z2": int(z2m.sum()), "Z3": int(z3m.sum()),
                       "total": int(len(positions))},
    "sequence_features": {
        "r_ev_cpg_oe": safe(pearsonr(ev_vals[gc_valid & ~np.isnan(cpg_arr)],
                                     cpg_arr[gc_valid & ~np.isnan(cpg_arr)])[0])
                        if chr1_seq else None,
        "zone_means_cpg_oe": {
            zone: safe(np.nanmean(cpg_arr[m & gc_valid]))
            for zone, m in [("Z1", z1m), ("Z2", z2m), ("Z3", z3m)]
        } if chr1_seq else {},
    },
    "exp1_methylation": {
        "r_ev_methylation":    safe(r_m)    if "r_m"    in dir() else None,
        "rho_ev_methylation":  safe(rho_m)  if "rho_m"  in dir() else None,
        "cohens_d_Z3_Z1":      safe(d_meth) if "d_meth" in dir() else None,
        "z1_mean_methylation": safe(np.mean(z1_m)) if "z1_m" in dir() else None,
        "z3_mean_methylation": safe(np.mean(z3_m)) if "z3_m" in dir() else None,
        "monotone_quartiles":  monotone if "monotone" in dir() else None,
    },
    "exp2_te": {
        "r_ev_te_fraction":  safe(r_te)  if "r_te"  in dir() else None,
        "rho_ev_te_fraction":safe(rho_te)if "rho_te" in dir() else None,
        "cohens_d_Z3_Z1":    safe(d_te)  if "d_te"   in dir() else None,
        "n_te_classes":      len(te_classes),
        "te_coverage_pct":   safe(te_cov_pct) if "te_cov_pct" in dir() else None,
        "te_class_r_with_ev": [
            {"class": cls, "r": safe(r), "rho": safe(rho)}
            for cls, r, rho in (cls_results if "cls_results" in dir() else [])
        ],
    },
    "exp3_cross_species": {
        "r_cpg_oe_vs_zone3pct":   safe(r_xs)   if "r_xs"  in dir() else None,
        "rho_cpg_oe_vs_zone3pct": safe(rho_xs) if "rho_xs" in dir() else None,
        "n_species":  len(xs_names),
        "species_data": {
            sp: {
                "cpg_oe":   sp_cpg_results[sp]["mean_cpg_oe"],
                "zone3_pct": float(xs_z3pct[xs_names.index(sp)])
                             if sp in xs_names else None,
                "clade":    xs_clades[xs_names.index(sp)]
                            if sp in xs_names else None,
            }
            for sp in sp_cpg_results
        },
    },
}

with open("methylation_te_results.json", "w") as f:
    json.dump(output, f, indent=2, default=safe)
print("  Saved: methylation_te_results.json\n")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# INTERPRETATION GUIDE
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
section("INTERPRETATION GUIDE")
print("""
  HOW TO READ RESULTS
  ───────────────────

  EXP 1 — Methylation:
    STRONG SUPPORT  Zone3 methylation >> Zone1, Cohen's d > 0.5
                    r(Ev, methylation) < −0.35
                    Monotonic quartile table (Q1→Q4 Ev strictly decreasing)
    WEAK / NULL     d < 0.2, r > −0.15 → methylation is NOT the mechanism
    RED FLAG        Zone1 has higher methylation than Zone3 → wrong direction,
                    fundamental error somewhere (wrong chr? wrong species?)

  EXP 2 — TE annotation:
    STRONG SUPPORT  OR(Zone3 in TE-rich) > 3, p < 1e-10
                    r(Ev, TE_frac) < −0.30
    WEAK / NULL     OR < 1.5 → TEs don't explain Zone3
    KEY FOLLOW-UP   Which TE class drives Zone3?
                    LTR retrotransposons (Gypsy/Copia) are the primary TEs
                    in plant centromeres/pericentromeres. If they dominate,
                    that closes the loop: LTR TEs → heterochromatin → CpG
                    suppression → low 4-mer skewness → Zone3.

  EXP 3 — Cross-species:
    STRONG SUPPORT  r(CpG_O/E, Zone3%) < −0.5 across ≥6 species
                    Monocots have lower CpG O/E than dicots
                    → monocot Zone3 shift = higher TE/methylation burden
    WEAK / NULL     r > −0.2 → inter-species Zone3 variation is NOT explained
                    by CpG suppression alone
    ALTERNATIVE     Zone3% correlates with TE% but not CpG O/E
                    → TE content is the proximate cause; methylation is secondary

  THE UNIFIED STORY (if all three pass):
    Ev_resid detects the euchromatin/heterochromatin boundary via its sensitivity
    to CpG dinucleotide frequency. Regions with high TE density are methylated
    → CpG suppression over evolutionary time → specific 4-mer depletion →
    lower-than-expected 4-mer skewness → Zone3.

    Gene-rich euchromatin is unmethylated → CpGs maintained → higher 4-mer
    skewness than GC content alone predicts → Zone1.

    The Zone3 mean (−2.83) is conserved at CV=2.6% across 1.5 Gyr because the
    fundamental TE-silencing methylation system is universal in eukaryotes,
    driving convergent CpG suppression to the same equilibrium regardless of
    TE family or genome size.
""")
