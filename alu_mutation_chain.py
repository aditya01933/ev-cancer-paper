#!/usr/bin/env python3
"""End-to-end: Alu density → mutation rate per window, chr1 BRCA."""
import numpy as np, json, gzip, os, sys
from scipy.stats import pearsonr, spearmanr

W = 5000
ZONES_MUT = os.path.expanduser("~/ai-project/ev-cancer-paper/raw_zones/human_chr1_ev_zones_WITH_MUTATIONS.json")
RMSK = os.path.expanduser("~/ai-project/ev-cancer-paper/raw_data/rmsk.txt.gz")

# === Load mutation-annotated zones ===
print("Loading mutation-annotated zones...")
with open(ZONES_MUT) as f:
    zones = json.load(f)
n_win = len(zones)
print(f"  {n_win} windows, keys: {list(zones[0].keys())}")

ev_resid = np.array([w.get('ev_resid', 0) for w in zones], dtype=float)
gc_vals = np.array([w.get('gc', 0) for w in zones], dtype=float)
win_starts = np.array([w['start'] for w in zones])
mut_total = np.array([w.get('mutation_count', 0) for w in zones], dtype=float)
print(f"  mut_total: nonzero={np.sum(mut_total>0)}, max={mut_total.max():.0f}, mean={mut_total.mean():.2f}")

# Zone labels
Z1_T, Z3_T = 0.382, -1.471
zone_lbl = np.where(ev_resid >= Z1_T, 'Z1', np.where(ev_resid <= Z3_T, 'Z3', 'Z2'))
print(f"  Zones: Z1={np.sum(zone_lbl=='Z1')}, Z2={np.sum(zone_lbl=='Z2')}, Z3={np.sum(zone_lbl=='Z3')}")

# === Build Alu density per window from RepeatMasker ===
print("\nParsing RepeatMasker for Alu density...")
alu_bp = np.zeros(n_win)
line_bp = np.zeros(n_win)
count = 0
with gzip.open(RMSK, 'rt') as f:
    for line in f:
        parts = line.split('\t')
        if len(parts) < 13 or parts[5] != 'chr1': continue
        start, end = int(parts[6]), int(parts[7])
        rep_class, rep_family = parts[11], parts[12]
        w_si = max(0, np.searchsorted(win_starts, start, side='right') - 1)
        w_ei = min(np.searchsorted(win_starts, end, side='right') - 1, n_win - 1)
        for wi in range(w_si, w_ei + 1):
            wl = win_starts[wi]
            wr = wl + W
            overlap = min(end, wr) - max(start, wl)
            if overlap <= 0: continue
            if rep_family == 'Alu':
                alu_bp[wi] += overlap
            elif rep_class == 'LINE':
                line_bp[wi] += overlap
        count += 1
        if count % 500000 == 0:
            print(f"  {count:,}...")
print(f"  Done: {count:,} records")

alu_frac = alu_bp / W

# === CORE TESTS ===
print("\n" + "=" * 60)
print("CORE TEST: Alu density → mutation rate")
print("=" * 60)
print(f"Windows with mutations: {np.sum(mut_total>0)}")

r, p = pearsonr(alu_bp, mut_total)
rs, ps = spearmanr(alu_bp, mut_total)
print(f"\nAlu_bp vs mutation_count (all {n_win} windows):")
print(f"  Pearson  r={r:+.4f}  p={p:.2e}")
print(f"  Spearman r={rs:+.4f}  p={ps:.2e}")

r2, p2 = pearsonr(ev_resid, mut_total)
print(f"\nEv_resid vs mutation_count:")
print(f"  Pearson  r={r2:+.4f}  p={p2:.2e}")

r3, p3 = pearsonr(line_bp, mut_total)
print(f"\nLINE_bp vs mutation_count:")
print(f"  Pearson  r={r3:+.4f}  p={p3:.2e}")

r4, p4 = pearsonr(gc_vals, mut_total)
print(f"\nGC vs mutation_count:")
print(f"  Pearson  r={r4:+.4f}  p={p4:.2e}")

# Per zone
print(f"\n{'Zone':<6} {'N':>6} {'Alu_mean':>10} {'Mut_mean':>10} {'MutRate/AluKb':>14}")
for z in ['Z1', 'Z2', 'Z3']:
    mask = zone_lbl == z
    n = mask.sum()
    am = alu_bp[mask].mean()
    mm = mut_total[mask].mean()
    rate = mm / max(am, 1) * 1000
    print(f"{z:<6} {n:6d} {am:10.1f} {mm:10.3f} {rate:14.4f}")

# === R² COMPARISON ===
print("\n" + "=" * 60)
print("R² COMPARISON: what predicts mutations best?")
print("=" * 60)
try:
    from sklearn.linear_model import LinearRegression
    target = mut_total

    r2_ev = LinearRegression().fit(ev_resid.reshape(-1,1), target).score(ev_resid.reshape(-1,1), target)
    r2_alu = LinearRegression().fit(alu_bp.reshape(-1,1), target).score(alu_bp.reshape(-1,1), target)
    r2_gc = LinearRegression().fit(gc_vals.reshape(-1,1), target).score(gc_vals.reshape(-1,1), target)
    r2_line = LinearRegression().fit(line_bp.reshape(-1,1), target).score(line_bp.reshape(-1,1), target)

    X_ev_alu = np.column_stack([ev_resid, alu_bp])
    r2_ev_alu = LinearRegression().fit(X_ev_alu, target).score(X_ev_alu, target)

    X_all = np.column_stack([ev_resid, alu_bp, line_bp, gc_vals])
    r2_all = LinearRegression().fit(X_all, target).score(X_all, target)

    print(f"  R² (Ev_resid alone):      {r2_ev:.4f}")
    print(f"  R² (Alu_bp alone):        {r2_alu:.4f}")
    print(f"  R² (GC alone):            {r2_gc:.4f}")
    print(f"  R² (LINE_bp alone):       {r2_line:.4f}")
    print(f"  R² (Ev_resid + Alu):      {r2_ev_alu:.4f}")
    print(f"  R² (Ev + Alu + LINE + GC):{r2_all:.4f}")
    print(f"\n  Alu adds beyond Ev:  {r2_ev_alu - r2_ev:.4f}")
    print(f"  Ev adds beyond Alu:  {r2_ev_alu - r2_alu:.4f}")
except ImportError:
    print("  sklearn not available")

print("\nDONE.")
