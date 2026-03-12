import json, numpy as np
from scipy.stats import fisher_exact

all_gc, all_mut = [], []

for chrom in [str(c) for c in range(1,23)] + ['X']:
    d = json.load(open(f'raw_zones/human_chr{chrom}_ev_zones_WITH_MUTATIONS.json'))
    all_gc.extend(w['gc'] for w in d)
    all_mut.extend(w['mutation_count'] for w in d)

all_gc = np.array(all_gc)
all_mut = np.array(all_mut)

assert len(all_gc) > 500000, f"Only {len(all_gc)} windows"
assert all_mut.sum() > 50000, f"Only {all_mut.sum()} mutations"

p_z1 = np.percentile(all_gc, 73.8)   # top 26.2%
p_z3 = np.percentile(all_gc, 24.4)   # bottom 24.4%

z1 = all_gc >= p_z1
z3 = all_gc <= p_z3

m_z3, w_z3 = all_mut[z3].sum(), z3.sum()
m_z1, w_z1 = all_mut[z1].sum(), z1.sum()

assert m_z3 > 0 and m_z1 > 0, "Zero mutations in a zone"

or_gc, p_gc = fisher_exact([[m_z3, w_z3-m_z3],[m_z1, w_z1-m_z1]], alternative='greater')

print(f"Windows: {len(all_gc):,}  |  Mutations: {all_mut.sum():,}")
print(f"GC thresholds: Z1>={p_z1:.4f}, Z3<={p_z3:.4f}")
print(f"GC-zones  OR={or_gc:.3f}  p={p_gc:.2e}  Z3_rate={m_z3/w_z3:.4f}  Z1_rate={m_z1/w_z1:.4f}")
print(f"Ev-zones  OR=1.682  (reference)")
print(f"Delta OR = {or_gc - 1.682:+.3f}")
