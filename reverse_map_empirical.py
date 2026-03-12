#!/usr/bin/env python3
"""Reverse-map Ev projection to identify driving 4-mers."""
import numpy as np
from scipy.stats import skew as skewness, pearsonr
from itertools import product
import json, sys, os

# === CANONICAL FORMULA ===
P = np.random.default_rng(42).standard_normal((256, 500)) / np.sqrt(256)
assert abs(P[0,0] - 0.01904482) < 1e-5, "WRONG RNG"
PRIMES = [p for p in range(2, 500) if all(p % i != 0 for i in range(2, p))][:95]
BASES = 'ACGT'
KMERS = [''.join(x) for x in product(BASES, repeat=4)]  # 256 4-mers
KI = {k: i for i, k in enumerate(KMERS)}
P_eff = P[:, PRIMES]  # 256 x 95 — the actual working projection

# === PART A: ANALYTICAL DECOMPOSITION (no genome needed) ===
print("=" * 60)
print("PART A: P-matrix weight decomposition")
print("=" * 60)

# For each 4-mer i, its 95-dim fingerprint is P_eff[i,:]
# Influence metrics:
#   1. L2 norm — total projection energy
#   2. Weight range (max-min) — spread across projections
#   3. Weight skewness — asymmetry of fingerprint (drives Ev skewness)
#   4. Mean abs weight — average influence per projection

results_a = []
for i, kmer in enumerate(KMERS):
    w = P_eff[i, :]
    has_cpg = 'CG' in kmer
    at_frac = sum(1 for c in kmer if c in 'AT') / 4
    results_a.append({
        'kmer': kmer, 'idx': i,
        'l2': np.linalg.norm(w),
        'wrange': float(w.max() - w.min()),
        'wskew': float(skewness(w)),
        'mean_abs': float(np.abs(w).mean()),
        'wmean': float(w.mean()),
        'wstd': float(w.std()),
        'has_cpg': has_cpg,
        'at_frac': at_frac
    })

# Rank by L2 norm (total influence)
results_a.sort(key=lambda x: x['l2'], reverse=True)

print(f"\nTop 20 4-mers by L2 projection weight:")
print(f"{'Rank':>4} {'4-mer':>6} {'L2':>7} {'wSkew':>7} {'wRange':>7} {'CpG':>4} {'AT%':>4}")
for r, d in enumerate(results_a[:20]):
    print(f"{r+1:4d} {d['kmer']:>6} {d['l2']:7.4f} {d['wskew']:+7.3f} {d['wrange']:7.4f} {'Y' if d['has_cpg'] else 'N':>4} {d['at_frac']:4.2f}")

print(f"\nBottom 10 (lowest influence):")
for r, d in enumerate(results_a[-10:]):
    print(f"{256-9+r:4d} {d['kmer']:>6} {d['l2']:7.4f} {d['wskew']:+7.3f} {d['wrange']:7.4f} {'Y' if d['has_cpg'] else 'N':>4} {d['at_frac']:4.2f}")

# CpG vs non-CpG group statistics
cpg_l2 = [d['l2'] for d in results_a if d['has_cpg']]
non_l2 = [d['l2'] for d in results_a if not d['has_cpg']]
print(f"\nCpG-containing 4-mers (n={len(cpg_l2)}): mean L2={np.mean(cpg_l2):.4f} ± {np.std(cpg_l2):.4f}")
print(f"Non-CpG 4-mers (n={len(non_l2)}):       mean L2={np.mean(non_l2):.4f} ± {np.std(non_l2):.4f}")

# AT-rich vs GC-rich
at_rich = [d['l2'] for d in results_a if d['at_frac'] >= 0.75]
gc_rich = [d['l2'] for d in results_a if d['at_frac'] <= 0.25]
mixed   = [d['l2'] for d in results_a if 0.25 < d['at_frac'] < 0.75]
print(f"AT-rich (≥3/4 AT, n={len(at_rich)}):     mean L2={np.mean(at_rich):.4f}")
print(f"GC-rich (≥3/4 GC, n={len(gc_rich)}):     mean L2={np.mean(gc_rich):.4f}")
print(f"Mixed (n={len(mixed)}):                   mean L2={np.mean(mixed):.4f}")

# === CRITICAL TEST: Does P_eff have systematic bias by composition? ===
# If random matrix is truly random, L2 norms should be ~equal for all 4-mers
# Any systematic difference = seed-specific artifact vs real structure
all_l2 = [d['l2'] for d in results_a]
print(f"\nL2 norm distribution: mean={np.mean(all_l2):.4f}, std={np.std(all_l2):.4f}, CV={np.std(all_l2)/np.mean(all_l2)*100:.1f}%")
print(f"Expected for random 95-dim Gaussian: L2 ≈ sqrt(95/256)*sqrt(95) = {np.sqrt(95/256)*np.sqrt(95):.4f}")

# === PART A2: Effective weight vector (how each 4-mer shifts mean of projections) ===
# w_eff[i] = mean(P_eff[i,:]) — if positive, increasing freq of kmer i raises mean projection
w_eff = P_eff.mean(axis=1)  # 256-vector
at_content = np.array([sum(1 for c in k if c in 'AT')/4 for k in KMERS])
r_at, p_at = pearsonr(w_eff, at_content)
gc_content = 1 - at_content
r_gc, p_gc = pearsonr(w_eff, gc_content)

cpg_flag = np.array([1 if 'CG' in k else 0 for k in KMERS])
r_cpg, p_cpg = pearsonr(w_eff, cpg_flag)

print(f"\nEffective weight vector correlations:")
print(f"  w_eff vs AT-content:  r={r_at:+.4f}, p={p_at:.2e}")
print(f"  w_eff vs GC-content:  r={r_gc:+.4f}, p={p_gc:.2e}")
print(f"  w_eff vs CpG-flag:    r={r_cpg:+.4f}, p={p_cpg:.2e}")

# Save Part A results
with open('reverse_map_analytical.json', 'w') as f:
    json.dump(results_a, f, indent=1)
print(f"\nPart A saved: reverse_map_analytical.json")

# === PART B: EMPIRICAL (needs chr1 FASTA + zone JSON) ===
print("\n" + "=" * 60)
print("PART B: Empirical correlation on chr1 windows")
print("=" * 60)

# Adjust these paths to your system
FASTA = os.path.expanduser("~/ai-project/ev-cancer-paper/raw_data/fasta/human_chr1.fa.gz")
ZONES = os.path.expanduser("~/ai-project/ev-cancer-paper/raw_zones/human_chr1_ev_zones.json")

if not os.path.exists(FASTA):
    print(f"FASTA not found: {FASTA}")
    print("Set FASTA path and rerun. Part A results above are still valid.")
    sys.exit(0)

# Read FASTA
import gzip
print("Reading chr1 FASTA (.gz)...")
seq = []
opener = gzip.open if FASTA.endswith('.gz') else open
with opener(FASTA, 'rt') as fh:
    for line in fh:
        if not line.startswith('>'): seq.append(line.strip().upper())
seq = ''.join(seq)
print(f"  Length: {len(seq):,} bp")

# Compute per-window 4-mer frequencies + Ev
W = 5000
n_windows = len(seq) // W
print(f"  Windows: {n_windows:,} ({W}bp)")

freqs = np.zeros((n_windows, 256), dtype=np.float64)
evs = np.zeros(n_windows)
gcs = np.zeros(n_windows)
valid = np.ones(n_windows, dtype=bool)

for wi in range(n_windows):
    s = seq[wi*W:(wi+1)*W]
    if 'N' in s and s.count('N') > W * 0.1:
        valid[wi] = False
        continue
    # 4-mer counts
    v = np.zeros(256)
    for i in range(len(s) - 3):
        k = s[i:i+4]
        if k in KI: v[KI[k]] += 1
    total = v.sum()
    if total < W * 0.5:
        valid[wi] = False
        continue
    f = v / total
    freqs[wi] = f
    gcs[wi] = sum(1 for c in s if c in 'GC') / len(s)
    ps = (f @ P)[PRIMES]
    evs[wi] = abs(skewness(ps)) * 6.07 + 0.10
    if (wi + 1) % 10000 == 0:
        print(f"  {wi+1}/{n_windows}...")

mask = valid
freqs, evs, gcs = freqs[mask], evs[mask], gcs[mask]
print(f"  Valid windows: {mask.sum():,}")

# Ev_resid
ev_resid = (evs - (-5.938 * gcs + 4.471)) / 0.456

# Correlate each 4-mer frequency with Ev_resid
print("\nComputing per-4-mer correlations with Ev_resid...")
results_b = []
for i, kmer in enumerate(KMERS):
    col = freqs[:, i]
    if col.std() < 1e-10:
        continue
    r, p = pearsonr(col, ev_resid)
    # Also partial: correlate with Ev_resid controlling for GC
    # (Ev_resid already residualized vs GC, so this IS the partial)
    results_b.append({
        'kmer': kmer, 'idx': i,
        'r_ev_resid': float(r), 'p_ev_resid': float(p),
        'mean_freq': float(col.mean()),
        'has_cpg': 'CG' in kmer,
        'at_frac': sum(1 for c in kmer if c in 'AT') / 4
    })

results_b.sort(key=lambda x: abs(x['r_ev_resid']), reverse=True)

print(f"\nTop 20 4-mers correlated with Ev_resid (GC-corrected):")
print(f"{'Rank':>4} {'4-mer':>6} {'r':>8} {'p':>10} {'MeanFreq':>9} {'CpG':>4} {'AT%':>4}")
for r, d in enumerate(results_b[:20]):
    print(f"{r+1:4d} {d['kmer']:>6} {d['r_ev_resid']:+8.4f} {d['p_ev_resid']:10.2e} {d['mean_freq']:9.6f} {'Y' if d['has_cpg'] else 'N':>4} {d['at_frac']:4.2f}")

print(f"\nBottom 10 (least correlated):")
for d in results_b[-10:]:
    print(f"     {d['kmer']:>6} {d['r_ev_resid']:+8.4f} {d['p_ev_resid']:10.2e} {'Y' if d['has_cpg'] else 'N':>4}")

# Summary: CpG vs non-CpG correlation magnitudes
cpg_r = [abs(d['r_ev_resid']) for d in results_b if d['has_cpg']]
non_r = [abs(d['r_ev_resid']) for d in results_b if not d['has_cpg']]
print(f"\nMean |r| with Ev_resid:")
print(f"  CpG-containing (n={len(cpg_r)}): {np.mean(cpg_r):.4f}")
print(f"  Non-CpG (n={len(non_r)}):        {np.mean(non_r):.4f}")
print(f"  Ratio: {np.mean(cpg_r)/np.mean(non_r):.2f}x")

# Top 5 positive + top 5 negative drivers
pos = [d for d in results_b if d['r_ev_resid'] > 0][:5]
neg = sorted([d for d in results_b if d['r_ev_resid'] < 0], key=lambda x: x['r_ev_resid'])[:5]
print(f"\nTop 5 POSITIVE drivers (high freq → high Ev_resid):")
for d in pos: print(f"  {d['kmer']} r={d['r_ev_resid']:+.4f} CpG={'Y' if d['has_cpg'] else 'N'}")
print(f"Top 5 NEGATIVE drivers (high freq → low Ev_resid):")
for d in neg: print(f"  {d['kmer']} r={d['r_ev_resid']:+.4f} CpG={'Y' if d['has_cpg'] else 'N'}")

# === PART B2: Can top-N 4-mers replace Ev? ===
# Build deterministic predictor from top 10 4-mers
top10_idx = [d['idx'] for d in results_b[:10]]
top10_signs = [np.sign(d['r_ev_resid']) for d in results_b[:10]]
# Simple weighted sum
det_score = freqs[:, top10_idx] @ np.array(top10_signs)
r_det, _ = pearsonr(det_score, ev_resid)
print(f"\nDeterministic 10-kmer score vs Ev_resid: r={r_det:.4f}")

# Compare to full Ev
from sklearn.metrics import r2_score  # noqa
try:
    from sklearn.linear_model import LinearRegression
    X_top = freqs[:, [d['idx'] for d in results_b[:20]]]
    lr = LinearRegression().fit(X_top, ev_resid)
    r2_20 = lr.score(X_top, ev_resid)
    print(f"Linear regression R² (top 20 4-mers → Ev_resid): {r2_20:.4f}")

    X_top5 = freqs[:, [d['idx'] for d in results_b[:5]]]
    lr5 = LinearRegression().fit(X_top5, ev_resid)
    r2_5 = lr5.score(X_top5, ev_resid)
    print(f"Linear regression R² (top 5 4-mers → Ev_resid):  {r2_5:.4f}")
except ImportError:
    print("sklearn not available — skip regression. pip install scikit-learn")

with open('reverse_map_empirical.json', 'w') as f:
    json.dump(results_b, f, indent=1)
print(f"\nPart B saved: reverse_map_empirical.json")
print("\nDONE. Send me both JSON files + terminal output for analysis.")
