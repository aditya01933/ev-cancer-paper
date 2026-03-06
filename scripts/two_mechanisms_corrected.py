"""
two_mechanisms_corrected.py — FIXED ChIP-seq two-mechanism test
Uses: BRCA GDC MAFs + fixed threshold zones + correct polarity
Z1 = high ev_resid >= +0.382 = euchromatin (GC-rich)
Z3 = low ev_resid <= -1.471 = heterochromatin (AT-rich)

Tests:
  Within Z3 (het): H3K9me3+ vs H3K9me3- mutation rate
  Within Z1 (eu):  H3K4me3+ vs H3K4me3- mutation rate

Available ChIP files:
  gm12878_h3k9me3_allchr.bed.gz  (genome-wide)
  gm12878_h3k4me3_chr1.bed.gz    (chr1 only)
  gm12878_h3k9me3_chr1.bed.gz    (chr1 only)

Run from: ~/ai-project/AncientKeyGen1/imp-research/tcga_cancer_experiment/
"""
import json, gzip, os, urllib.request
import numpy as np
from scipy.stats import mannwhitneyu
from collections import defaultdict

os.chdir(os.path.expanduser("~/ai-project/AncientKeyGen1/imp-research/tcga_cancer_experiment"))

# ── Verified constants ──
Z3_T, Z1_T = -1.471, 0.382
P = np.random.default_rng(42).standard_normal((256,500)) / np.sqrt(256)
assert abs(P[0,0] - 0.01904482) < 1e-5, "FATAL: Wrong RNG"
WINDOW = 5000

# ══════════════════════════════════════════════════════════════════
# STEP 1: Load zone files (base, NOT _WITH_MUTATIONS)
# ══════════════════════════════════════════════════════════════════
print("═══ STEP 1: Load zones (fixed thresholds) ═══")
zones = {}  # {(chrom_str, start): {'er': float, 'zone': str, 'muts': 0, 'k9': 0, 'k4': 0}}
zone_wins = defaultdict(int)

for c in list(range(1,23)) + ['X']:
    fname = f"human_chr{c}_ev_zones.json"
    if not os.path.exists(fname): continue
    with open(fname) as f: data = json.load(f)
    assert 'ev_resid' in data[0], f"{fname} missing ev_resid"
    for w in data:
        er = w['ev_resid']
        z = 'Z3' if er <= Z3_T else ('Z1' if er >= Z1_T else 'Z2')
        zones[(str(c), w['start'])] = {'er': er, 'zone': z, 'muts': 0, 'k9': 0, 'k4': 0}
        zone_wins[z] += 1

print(f"  Total windows: {len(zones)}")
print(f"  Z1(eu)={zone_wins['Z1']}  Z2={zone_wins['Z2']}  Z3(het)={zone_wins['Z3']}")
assert zone_wins['Z1'] > 100000 and zone_wins['Z3'] > 100000, "Zone counts wrong"

# ══════════════════════════════════════════════════════════════════
# STEP 2: Load ChIP-seq peaks → mark windows
# ══════════════════════════════════════════════════════════════════
print("\n═══ STEP 2: Load ChIP-seq peaks ═══")

def load_bed_peaks(bed_path):
    """Load BED file, return {(chrom, start_5kb): 1} for overlapping windows."""
    hits = {}
    opener = gzip.open if bed_path.endswith('.gz') else open
    with opener(bed_path, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'): continue
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            chrom = parts[0].replace('chr','')
            try:
                pstart, pend = int(parts[1]), int(parts[2])
            except ValueError: continue
            # Mark all 5kb windows overlapping this peak
            w_start = (pstart // WINDOW) * WINDOW
            w_end = (pend // WINDOW) * WINDOW
            for ws in range(w_start, w_end + WINDOW, WINDOW):
                hits[(chrom, ws)] = 1
    return hits

# H3K9me3 — genome-wide available
k9_file = "gm12878_h3k9me3_allchr.bed.gz"
if not os.path.exists(k9_file):
    k9_file = "gm12878_h3k9me3_chr1.bed.gz"
    print(f"  WARNING: using chr1-only H3K9me3")
k9_peaks = load_bed_peaks(k9_file)
print(f"  H3K9me3 peaks loaded from {k9_file}: {len(k9_peaks)} window overlaps")

# H3K4me3 — chr1 only available
k4_file = "gm12878_h3k4me3_chr1.bed.gz"
if os.path.exists("gm12878_h3k4me3_allchr.bed.gz"):
    k4_file = "gm12878_h3k4me3_allchr.bed.gz"
k4_peaks = load_bed_peaks(k4_file)
print(f"  H3K4me3 peaks loaded from {k4_file}: {len(k4_peaks)} window overlaps")

# Mark windows
n_k9 = n_k4 = 0
for key, w in zones.items():
    chrom, start = key
    if (chrom, start) in k9_peaks:
        w['k9'] = 1; n_k9 += 1
    if (chrom, start) in k4_peaks:
        w['k4'] = 1; n_k4 += 1

print(f"  Windows with H3K9me3: {n_k9}")
print(f"  Windows with H3K4me3: {n_k4}")

# SANITY: H3K9me3 should be more in Z3(het), H3K4me3 more in Z1(eu)
z3_k9 = sum(1 for k,w in zones.items() if w['zone']=='Z3' and w['k9'])
z1_k9 = sum(1 for k,w in zones.items() if w['zone']=='Z1' and w['k9'])
z3_k4 = sum(1 for k,w in zones.items() if w['zone']=='Z3' and w['k4'])
z1_k4 = sum(1 for k,w in zones.items() if w['zone']=='Z1' and w['k4'])
print(f"\n  SANITY — H3K9me3: Z3(het)={z3_k9}, Z1(eu)={z1_k9}")
print(f"  SANITY — H3K4me3: Z3(het)={z3_k4}, Z1(eu)={z1_k4}")

# ══════════════════════════════════════════════════════════════════
# STEP 3: Count BRCA mutations per window (from GDC, not MC3)
# ══════════════════════════════════════════════════════════════════
print("\n═══ STEP 3: Count BRCA mutations from GDC ═══")
hits_maf = json.load(open("maf_ids_TCGA_BRCA.json"))
print(f"  BRCA MAF files: {len(hits_maf)}")

n_files = n_mapped = 0
for hit in hits_maf:
    try:
        url = f"https://api.gdc.cancer.gov/data/{hit['file_id']}"
        req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as r: raw = r.read()
        try: content = gzip.decompress(raw).decode('utf-8',errors='replace')
        except: content = raw.decode('utf-8',errors='replace')
        lines = content.split('\n')
        header = None; ds = 0
        for i, line in enumerate(lines):
            if line.startswith('#'): continue
            header = line.split('\t'); ds = i+1; break
        if not header: continue
        chr_c = next((i for i,h in enumerate(header) if h=='Chromosome'), None)
        pos_c = next((i for i,h in enumerate(header) if h=='Start_Position'), None)
        if chr_c is None or pos_c is None: continue
        for line in lines[ds:]:
            if not line.strip(): continue
            p = line.split('\t')
            if len(p) <= max(chr_c, pos_c): continue
            chrom = p[chr_c].replace('chr','')
            start = (int(p[pos_c]) // WINDOW) * WINDOW
            key = (chrom, start)
            if key in zones:
                zones[key]['muts'] += 1
                n_mapped += 1
        n_files += 1
        if n_files % 100 == 0:
            print(f"    {n_files}/{len(hits_maf)} files, {n_mapped} mutations mapped")
    except: continue

print(f"  Done: {n_files} files, {n_mapped} mutations")

# SANITY: total should match verified result
z1_m = sum(w['muts'] for w in zones.values() if w['zone']=='Z1')
z3_m = sum(w['muts'] for w in zones.values() if w['zone']=='Z3')
print(f"  Z1(eu) muts={z1_m}, Z3(het) muts={z3_m}")
print(f"  Z3/Z1 rate = {(z3_m/zone_wins['Z3'])/(z1_m/zone_wins['Z1']):.4f} (expect ~1.55)")
assert abs(z1_m - 18413) < 200, f"Z1 muts mismatch: {z1_m}"
assert abs(z3_m - 26605) < 200, f"Z3 muts mismatch: {z3_m}"
print(f"  ✓ Matches verified OR=1.682 result")

# ══════════════════════════════════════════════════════════════════
# STEP 4: Two-mechanism test
# ══════════════════════════════════════════════════════════════════
print("\n═══ STEP 4: Two-mechanism analysis ═══")

def mechanism_test(zone_label, mark_field, mark_name):
    """Within a zone, compare mutation rate in mark+ vs mark- windows."""
    plus  = [w['muts'] for w in zones.values() if w['zone']==zone_label and w[mark_field]==1]
    minus = [w['muts'] for w in zones.values() if w['zone']==zone_label and w[mark_field]==0]
    if len(plus) < 10 or len(minus) < 10:
        print(f"  SKIP: {mark_name} in {zone_label} — too few windows ({len(plus)} plus, {len(minus)} minus)")
        return None
    mean_p, mean_m = np.mean(plus), np.mean(minus)
    ratio = mean_p / mean_m if mean_m > 0 else float('inf')
    U, p = mannwhitneyu(plus, minus, alternative='two-sided')
    print(f"\n  {zone_label} + {mark_name}:")
    print(f"    {mark_name}+ : n={len(plus):,}  mean_muts={mean_p:.4f}")
    print(f"    {mark_name}- : n={len(minus):,}  mean_muts={mean_m:.4f}")
    print(f"    Ratio ({mark_name}+ / {mark_name}-) = {ratio:.3f}")
    print(f"    Mann-Whitney p = {p:.2e}")
    if ratio > 1 and p < 0.05:
        print(f"    → {mark_name}+ windows have MORE mutations ✓")
    elif ratio < 1 and p < 0.05:
        print(f"    → {mark_name}+ windows have FEWER mutations ✗")
    else:
        print(f"    → Not significant")
    return {'n_plus': len(plus), 'n_minus': len(minus),
            'mean_plus': round(mean_p,4), 'mean_minus': round(mean_m,4),
            'ratio': round(ratio,4), 'p': float(p)}

print("\n── Mechanism 1: Heterochromatin repair failure ──")
print("  Hypothesis: Z3(het) + H3K9me3+ has MORE mutations than Z3(het) + H3K9me3-")
m1 = mechanism_test('Z3', 'k9', 'H3K9me3')

print("\n── Mechanism 2: Euchromatin transcription-coupled ──")
print("  Hypothesis: Z1(eu) + H3K4me3+ has MORE mutations than Z1(eu) + H3K4me3-")
print("  (NOTE: H3K4me3 is chr1-only unless allchr file found)")
m2 = mechanism_test('Z1', 'k4', 'H3K4me3')

# Also test cross-zone ChIP effects
print("\n── Additional: Z3(het) + H3K4me3 ──")
m3 = mechanism_test('Z3', 'k4', 'H3K4me3')

print("\n── Additional: Z1(eu) + H3K9me3 ──")
m4 = mechanism_test('Z1', 'k9', 'H3K9me3')

# ══════════════════════════════════════════════════════════════════
# STEP 5: Summary
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("SUMMARY — Two-Mechanism Test (CORRECTED)")
print("="*60)
print(f"Data: BRCA GDC ({n_files} patients), {n_mapped} mutations")
print(f"Zones: Z3=het (ev_resid≤{Z3_T}), Z1=eu (ev_resid≥{Z1_T})")
print(f"ChIP: GM12878 H3K9me3 ({k9_file}), H3K4me3 ({k4_file})")
print()
if m1:
    tag = "CONFIRMED" if m1['ratio']>1 and m1['p']<0.05 else "NOT CONFIRMED"
    print(f"Mech 1 (het repair failure): ratio={m1['ratio']}, p={m1['p']:.2e} → {tag}")
if m2:
    tag = "CONFIRMED" if m2['ratio']>1 and m2['p']<0.05 else "NOT CONFIRMED"
    print(f"Mech 2 (eu transcription):   ratio={m2['ratio']}, p={m2['p']:.2e} → {tag}")

results = {'mechanism_1_het_repair': m1, 'mechanism_2_eu_transcription': m2,
           'z3_h3k4me3': m3, 'z1_h3k9me3': m4,
           'data': {'n_files': n_files, 'n_mapped': n_mapped, 'z1_muts': z1_m, 'z3_muts': z3_m}}
with open('two_mechanisms_corrected.json','w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved: two_mechanisms_corrected.json")
