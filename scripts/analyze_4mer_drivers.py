"""
Which 4-mers drive Ev_resid after GC removal?
Tests nucleosome-positioning hypothesis (AA/TT at 10bp periodicity).
Uses chr1 only. Runtime ~5 min.
"""
import gzip, json, numpy as np, sys
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict

sys.path.insert(0,'scripts')
from ev_formula import KI, Z1_T, Z3_T

FASTA = 'raw_data/fasta/human_chr1.fa.gz'
ZONES = 'raw_zones/human_chr1_ev_zones_WITH_MUTATIONS.json'
WS    = 5000
BASES = 'ACGT'
ALL_4MERS = [a+b+c+d for a in BASES for b in BASES for c in BASES for d in BASES]

# ── Load zones ──
zones = json.load(open(ZONES))
assert all('gc' in w and 'ev_resid' in w for w in zones), "Missing gc or ev_resid"
starts = {w['start']: w for w in zones}
print(f"Loaded {len(zones)} windows")

# ── Load FASTA ──
with gzip.open(FASTA,'rt') as f:
    seq = ''.join(l.strip() for l in f if not l.startswith('>'))
assert len(seq) > 1_000_000
print(f"Loaded {len(seq):,} bp")

# ── Compute 4-mer frequencies per window ──
print("Computing 4-mer frequencies...", flush=True)
ev_vals, gc_vals, kmer_mat = [], [], []

for w in zones:
    seg = seq[w['start']:w['start']+WS].upper()
    seg_clean = seg.replace('N','')
    if len(seg_clean) < 1000: continue

    v = np.zeros(256)
    for i in range(len(seg)-3):
        k = seg[i:i+4]
        if k in KI: v[KI[k]] += 1
    if v.sum() == 0: continue

    f = v / v.sum()
    ev_vals.append(w['ev_resid'])
    gc_vals.append(w['gc'])
    kmer_mat.append(f)

ev_arr = np.array(ev_vals)
gc_arr = np.array(gc_vals)
K      = np.array(kmer_mat)  # shape: (n_windows, 256)

assert K.shape[1] == 256
assert len(ev_arr) > 10000, f"Too few windows: {len(ev_arr)}"
print(f"Matrix: {K.shape}  ev_mean={ev_arr.mean():.3f}")

# ── Partial correlation: each 4-mer vs ev_resid controlling for GC ──
# Residualize both kmer_freq and ev_resid on GC, then correlate
from numpy.linalg import lstsq

def residualize(y, X):
    """Remove linear effect of X from y."""
    X2 = np.column_stack([np.ones(len(X)), X])
    coef, *_ = lstsq(X2, y, rcond=None)
    return y - X2 @ coef

ev_resid_on_gc = residualize(ev_arr, gc_arr)

print("\nTop 4-mers correlated with Ev_resid AFTER removing GC effect:")
print(f"{'4-mer':>6}  {'r(raw,ev)':>10}  {'r(gc-removed,ev)':>18}  {'Z1_freq':>8}  {'Z3_freq':>8}  {'Z3/Z1':>6}")

results = []
z1_mask = ev_arr >= Z1_T
z3_mask = ev_arr <= Z3_T

for i, kmer in enumerate(ALL_4MERS):
    freq = K[:, i]
    freq_resid = residualize(freq, gc_arr)

    r_raw,  _ = pearsonr(freq, ev_arr)
    r_resid,_ = pearsonr(freq_resid, ev_resid_on_gc)

    z1_f = K[z1_mask, i].mean()
    z3_f = K[z3_mask, i].mean()
    ratio = z3_f / z1_f if z1_f > 0 else np.nan

    results.append((kmer, r_raw, r_resid, z1_f, z3_f, ratio))

# Sort by abs partial correlation
results.sort(key=lambda x: abs(x[2]), reverse=True)

print("\n--- TOP 20 4-mers driving Ev beyond GC ---")
for kmer, r_raw, r_resid, z1_f, z3_f, ratio in results[:20]:
    print(f"{kmer:>6}  {r_raw:>+10.4f}  {r_resid:>+18.4f}  {z1_f:>8.5f}  {z3_f:>8.5f}  {ratio:>6.3f}")

# ── Nucleosome positioning test: AA/TT dinucleotide periodicity ──
print("\n--- Nucleosome positioning 4-mers (AA/TT rich) ---")
nuc_kmers = [k for k in ALL_4MERS if k.count('A')+k.count('T') >= 3]
nuc_results = [(k,r,p,z1,z3,rat) for k,r,p,z1,z3,rat in results if k in nuc_kmers]
nuc_results.sort(key=lambda x: abs(x[2]), reverse=True)
for kmer, r_raw, r_resid, z1_f, z3_f, ratio in nuc_results[:10]:
    print(f"{kmer:>6}  {r_raw:>+10.4f}  {r_resid:>+18.4f}  {z1_f:>8.5f}  {z3_f:>8.5f}  {ratio:>6.3f}")

# ── Save results ──
out = {
    'n_windows': len(ev_arr),
    'top20_partial_corr': [
        {'kmer':k,'r_raw':round(r,4),'r_partial':round(p,4),
         'z1_freq':round(z1,6),'z3_freq':round(z3,6),'z3_z1_ratio':round(rat,4)}
        for k,r,p,z1,z3,rat in results[:20]
    ],
    'top10_AT_rich': [
        {'kmer':k,'r_raw':round(r,4),'r_partial':round(p,4),
         'z1_freq':round(z1,6),'z3_freq':round(z3,6),'z3_z1_ratio':round(rat,4)}
        for k,r,p,z1,z3,rat in nuc_results[:10]
    ]
}
json.dump(out, open('data/ev_4mer_drivers.json','w'), indent=2)
print(f"\nSaved → data/ev_4mer_drivers.json")
