# scripts/test_at_skew.py
import gzip, json, numpy as np
from scipy.stats import pearsonr, spearmanr

FASTA = 'raw_data/fasta/human_chr1.fa.gz'
ZONES = 'raw_zones/human_chr1_ev_zones_WITH_MUTATIONS.json'
WS    = 5000

zones = json.load(open(ZONES))
with gzip.open(FASTA,'rt') as f:
    seq = ''.join(l.strip() for l in f if not l.startswith('>'))

ev_vals, gc_vals, at_skew, cg_skew = [], [], [], []

for w in zones:
    seg = seq[w['start']:w['start']+WS].upper().replace('N','')
    if len(seg) < 1000: continue
    A,T,G,C = seg.count('A'),seg.count('T'),seg.count('G'),seg.count('C')
    total = A+T+G+C
    if total == 0: continue
    ev_vals.append(w['ev_resid'])
    gc_vals.append(w['gc'])
    at_skew.append((A-T)/max(A+T,1))   # strand asymmetry: A vs T
    cg_skew.append((C-G)/max(C+G,1))   # strand asymmetry: C vs G

ev_arr = np.array(ev_vals)
gc_arr = np.array(gc_vals)
at_arr = np.array(at_skew)
cg_arr = np.array(cg_skew)

# Residualize ev and skews on GC
from numpy.linalg import lstsq
def resid(y, x):
    X = np.column_stack([np.ones(len(x)), x])
    return y - X @ lstsq(X, y, rcond=None)[0]

ev_r  = resid(ev_arr, gc_arr)
at_r  = resid(at_arr, gc_arr)
cg_r  = resid(cg_arr, gc_arr)

r_at_raw,  _ = pearsonr(at_arr, ev_arr)
r_at_resid,_ = pearsonr(at_r,   ev_r)
r_cg_raw,  _ = pearsonr(cg_arr, ev_arr)
r_cg_resid,_ = pearsonr(cg_r,   ev_r)

print(f"AT-skew vs Ev:         r_raw={r_at_raw:+.4f}  r_partial={r_at_resid:+.4f}")
print(f"CG-skew vs Ev:         r_raw={r_cg_raw:+.4f}  r_partial={r_cg_resid:+.4f}")

# By zone
z1 = ev_arr >= 0.382
z3 = ev_arr <= -1.471
print(f"\nAT-skew: Z1={at_arr[z1].mean():+.4f}  Z3={at_arr[z3].mean():+.4f}  diff={at_arr[z1].mean()-at_arr[z3].mean():+.4f}")
print(f"CG-skew: Z1={cg_arr[z1].mean():+.4f}  Z3={cg_arr[z3].mean():+.4f}  diff={cg_arr[z1].mean()-cg_arr[z3].mean():+.4f}")

# Is AT-skew itself predictive of mutation rate?
mut = np.array([w['mutation_count'] for w in zones if len(seq[w['start']:w['start']+WS].upper().replace('N',''))>=1000])
r_mut_at,_ = pearsonr(at_arr[:len(mut)], mut)
r_mut_ev,_ = pearsonr(ev_arr[:len(mut)], mut)
print(f"\nMutation correlation: r(AT-skew,mut)={r_mut_at:+.4f}  r(Ev,mut)={r_mut_ev:+.4f}")

json.dump({
    'r_at_raw': round(r_at_raw,4), 'r_at_partial': round(r_at_resid,4),
    'r_cg_raw': round(r_cg_raw,4), 'r_cg_partial': round(r_cg_resid,4),
    'z1_at_skew': round(float(at_arr[z1].mean()),4),
    'z3_at_skew': round(float(at_arr[z3].mean()),4),
    'r_mut_at': round(r_mut_at,4), 'r_mut_ev': round(r_mut_ev,4)
}, open('data/at_skew_analysis.json','w'), indent=2)
print("\nSaved → data/at_skew_analysis.json")
