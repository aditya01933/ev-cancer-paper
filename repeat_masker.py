#!/usr/bin/env python3
"""Correlate Ev_resid with RepeatMasker TE families per 5kb window on chr1."""
import numpy as np, json, gzip, sys, os
from scipy.stats import pearsonr
from collections import defaultdict

W = 5000
ZONES = os.path.expanduser("~/ai-project/ev-cancer-paper/raw_zones/human_chr1_ev_zones.json")
RMSK = os.path.expanduser("~/ai-project/ev-cancer-paper/raw_data/rmsk.txt.gz")

if not os.path.exists(RMSK):
    print(f"RepeatMasker not found at {RMSK}")
    print("Download: curl -o ~/ai-project/ev-cancer-paper/raw_data/rmsk.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz")
    sys.exit(1)

# === Load zone data ===
print("Loading zone JSON...")
with open(ZONES) as f:
    zones = json.load(f)

n_win = len(zones)
ev_resid = np.array([w['ev_resid'] for w in zones])
win_starts = np.array([w['start'] for w in zones])

# Assign zones from ev_resid using canonical thresholds
Z1_T, Z3_T = 0.382, -1.471
zone_lbl = np.where(ev_resid >= Z1_T, 'Z1', np.where(ev_resid <= Z3_T, 'Z3', 'Z2'))
print(f"  {n_win} windows")
print(f"  Zone counts: Z1={np.sum(zone_lbl=='Z1')}, Z2={np.sum(zone_lbl=='Z2')}, Z3={np.sum(zone_lbl=='Z3')}")

valid = np.ones(n_win, dtype=bool)

# === Parse RepeatMasker for chr1 ===
print("Parsing RepeatMasker (chr1 only)...")
te_classes = ['LINE', 'SINE', 'LTR', 'DNA', 'Satellite', 'Simple_repeat', 'Low_complexity']
win_te_bp = {c: np.zeros(n_win) for c in te_classes}
win_line_sub = defaultdict(lambda: np.zeros(n_win))
win_sine_sub = defaultdict(lambda: np.zeros(n_win))
win_sat_sub  = defaultdict(lambda: np.zeros(n_win))

count = 0
with gzip.open(RMSK, 'rt') as f:
    for line in f:
        parts = line.split('\t')
        if len(parts) < 13: continue
        if parts[5] != 'chr1': continue
        start, end = int(parts[6]), int(parts[7])
        rep_class, rep_family, rep_name = parts[11], parts[12], parts[10]

        w_si = max(0, np.searchsorted(win_starts, start, side='right') - 1)
        w_ei = min(np.searchsorted(win_starts, end, side='right') - 1, n_win - 1)
        for wi in range(w_si, w_ei + 1):
            wl = win_starts[wi]
            wr = wl + W
            overlap = min(end, wr) - max(start, wl)
            if overlap <= 0: continue
            if rep_class in win_te_bp:
                win_te_bp[rep_class][wi] += overlap
            if rep_class == 'LINE':
                win_line_sub[rep_family][wi] += overlap
            elif rep_class == 'SINE':
                win_sine_sub[rep_family][wi] += overlap
            elif rep_class == 'Satellite':
                win_sat_sub[rep_family][wi] += overlap
        count += 1
        if count % 500000 == 0:
            print(f"  {count:,} records...")

print(f"  Total chr1 records: {count:,}")

# === Correlations ===
print("\n" + "=" * 60)
print("TE CLASS vs Ev_resid (GC-corrected)")
print("=" * 60)
print(f"{'Class':<20} {'r':>8} {'p':>12} {'MeanBP/win':>12} {'Z1_mean':>10} {'Z3_mean':>10}")
for cls in te_classes:
    arr = win_te_bp[cls]
    if arr.std() < 1e-10: continue
    r, p = pearsonr(arr, ev_resid)
    z1m = arr[zone_lbl == 'Z1'].mean() if np.sum(zone_lbl == 'Z1') > 0 else 0
    z3m = arr[zone_lbl == 'Z3'].mean() if np.sum(zone_lbl == 'Z3') > 0 else 0
    print(f"{cls:<20} {r:+8.4f} {p:12.2e} {arr.mean():12.1f} {z1m:10.1f} {z3m:10.1f}")

print("\n" + "=" * 60)
print("LINE subfamilies vs Ev_resid")
print("=" * 60)
for fam, arr in sorted(win_line_sub.items(), key=lambda x: -x[1].mean()):
    if arr.std() < 1e-10 or arr.mean() < 10: continue
    r, p = pearsonr(arr, ev_resid)
    print(f"  {fam:<15} r={r:+.4f}  p={p:.2e}  mean={arr.mean():.0f}bp/win")

print("\n" + "=" * 60)
print("SINE subfamilies vs Ev_resid")
print("=" * 60)
for fam, arr in sorted(win_sine_sub.items(), key=lambda x: -x[1].mean()):
    if arr.std() < 1e-10 or arr.mean() < 10: continue
    r, p = pearsonr(arr, ev_resid)
    print(f"  {fam:<15} r={r:+.4f}  p={p:.2e}  mean={arr.mean():.0f}bp/win")

print("\n" + "=" * 60)
print("Satellite subfamilies vs Ev_resid")
print("=" * 60)
for fam, arr in sorted(win_sat_sub.items(), key=lambda x: -x[1].mean()):
    if arr.std() < 1e-10 or arr.mean() < 5: continue
    r, p = pearsonr(arr, ev_resid)
    print(f"  {fam:<15} r={r:+.4f}  p={p:.2e}  mean={arr.mean():.0f}bp/win")

# === Zone summary ===
print("\n" + "=" * 60)
print("ZONE SUMMARY — mean TE bp per 5kb window")
print("=" * 60)
print(f"{'TE Class':<20} {'Z1':>10} {'Z2':>10} {'Z3':>10} {'Z3/Z1':>8}")
for cls in te_classes:
    arr = win_te_bp[cls]
    z1 = arr[zone_lbl == 'Z1'].mean() if np.sum(zone_lbl == 'Z1') > 0 else 0
    z2 = arr[zone_lbl == 'Z2'].mean() if np.sum(zone_lbl == 'Z2') > 0 else 0
    z3 = arr[zone_lbl == 'Z3'].mean() if np.sum(zone_lbl == 'Z3') > 0 else 0
    ratio = z3 / z1 if z1 > 0 else float('inf')
    print(f"{cls:<20} {z1:10.1f} {z2:10.1f} {z3:10.1f} {ratio:8.2f}")

print("\nDONE.")
