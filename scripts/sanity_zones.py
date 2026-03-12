# scripts/sanity_zones.py
import json, numpy as np

chroms = [str(c) for c in range(1,23)] + ['X']
total_w, total_m = 0, 0

for chrom in chroms:
    d = json.load(open(f'raw_zones/human_chr{chrom}_ev_zones_WITH_MUTATIONS.json'))
    keys = set(d[0].keys())
    ev   = [w['ev_resid'] for w in d]
    mut  = [w['mutation_count'] for w in d]
    gc   = [w['gc'] for w in d] if 'gc' in d[0] else None

    n_z1 = sum(1 for e in ev if e >= 0.382)
    n_z3 = sum(1 for e in ev if e <= -1.471)
    m_sum = sum(mut)
    neg   = sum(1 for m in mut if m < 0)

    total_w += len(d); total_m += m_sum
    gc_str = f"gc=[{min(gc):.3f},{max(gc):.3f}]" if gc else "gc=MISSING"
    print(f"chr{chrom:>2}: {len(d):>6}w  mut={m_sum:>7}  Z1={n_z1:>5}  Z3={n_z3:>5}  neg_mut={neg}  {gc_str}  keys={keys}")

print(f"\nTOTAL: {total_w:,} windows  {total_m:,} mutations")
assert total_w > 500000, "FAIL: too few windows"
assert total_m > 50000,  "FAIL: too few mutations"
print("PASS: counts look sane")
