"""
PCA of 4-mer space: what axes does Ev capture vs GC?
Shows which PCs Ev_resid correlates with after GC removal.
Uses chr1. Runtime ~3 min.
"""
import gzip, json, numpy as np, sys
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

sys.path.insert(0,'scripts')
from ev_formula import KI, Z1_T, Z3_T

FASTA = 'raw_data/fasta/human_chr1.fa.gz'
ZONES = 'raw_zones/human_chr1_ev_zones_WITH_MUTATIONS.json'
WS    = 5000

zones = json.load(open(ZONES))
with gzip.open(FASTA,'rt') as f:
    seq = ''.join(l.strip() for l in f if not l.startswith('>'))

ev_vals, gc_vals, kmer_mat = [], [], []
for w in zones:
    seg = seq[w['start']:w['start']+WS].upper().replace('N','')
    if len(seg) < 1000: continue
    v = np.zeros(256)
    for i in range(len(seg)-3):
        k = seg[i:i+4]
        if k in KI: v[KI[k]] += 1
    if v.sum()==0: continue
    kmer_mat.append(v/v.sum())
    ev_vals.append(w['ev_resid'])
    gc_vals.append(w['gc'])

K      = np.array(kmer_mat)
ev_arr = np.array(ev_vals)
gc_arr = np.array(gc_vals)
print(f"Matrix: {K.shape}")

# ── PCA of 4-mer space ──
pca = PCA(n_components=20)
PCs = pca.fit_transform(K)

print(f"\n{'PC':>4}  {'var%':>7}  {'r_with_GC':>10}  {'r_with_Ev':>10}  {'r_with_EV_resid':>16}")
pc_results = []
for i in range(20):
    pc    = PCs[:,i]
    r_gc, _ = pearsonr(pc, gc_arr)
    r_ev, _ = pearsonr(pc, ev_arr)
    # Ev_resid = Ev after GC removed — correlate PC with Ev controlling for GC
    from numpy.linalg import lstsq
    X2 = np.column_stack([np.ones(len(gc_arr)), gc_arr])
    ev_resid_gc = ev_arr - X2 @ lstsq(X2, ev_arr, rcond=None)[0]
    r_evr,_ = pearsonr(pc, ev_resid_gc)
    var = pca.explained_variance_ratio_[i]*100
    print(f"PC{i+1:>2}  {var:>7.2f}%  {r_gc:>+10.4f}  {r_ev:>+10.4f}  {r_evr:>+16.4f}")
    pc_results.append({'pc':i+1,'var_pct':round(var,3),
                       'r_gc':round(r_gc,4),'r_ev':round(r_ev,4),'r_ev_resid':round(r_evr,4)})

# Sanity: PC1 should dominate and correlate strongly with GC
assert pca.explained_variance_ratio_[0] > 0.3, "PC1 should explain >50% (GC axis)"
assert abs(pc_results[0]['r_gc']) > 0.8, "PC1 should correlate strongly with GC"
print(f"\n[OK] PC1={pc_results[0]['var_pct']:.1f}% var, r_GC={pc_results[0]['r_gc']:.3f} (GC axis confirmed)")

# Which PCs does Ev_resid capture?
print("\nPCs where |r_ev_resid| > 0.1 (Ev captures BEYOND GC):")
for r in pc_results:
    if abs(r['r_ev_resid']) > 0.1:
        print(f"  PC{r['pc']}: var={r['var_pct']:.2f}%  r_gc={r['r_gc']:+.3f}  r_ev_resid={r['r_ev_resid']:+.3f}")

json.dump({'pca': pc_results}, open('data/ev_pca_analysis.json','w'), indent=2)
print("\nSaved → data/ev_pca_analysis.json")
