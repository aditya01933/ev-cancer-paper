#!/usr/bin/env python3
"""tsg_onc_validation.py — 6-experiment TSG/oncogene zone segregation validation"""
import json, numpy as np
from scipy.stats import mannwhitneyu, fisher_exact, spearmanr, chi2_contingency

GENE_FILE = "/Users/aditya/ai-project/AncientKeyGen1/imp-research/tcga_cancer_experiment/drug_targets_zone3.json"
OUT = "/Users/aditya/ai-project/AncientKeyGen1/imp-research/tcga_cancer_experiment/tsg_onc_validation.json"

with open(GENE_FILE) as f:
    genes = json.load(f)

# Sanity
assert all('mean_ev' in g for g in genes)
assert all('role' in g for g in genes)
roles_seen = set(g['role'] for g in genes)
print(f"Roles in data: {roles_seen}")  # show actual values before filtering

tsgs = [g for g in genes if g['role'].upper() == 'TSG']
oncs = [g for g in genes if g['role'].upper() == 'ONCOGENE']
print(f"TSGs={len(tsgs)}, Oncogenes={len(oncs)} (of {len(genes)} total)")
assert len(tsgs) > 0 and len(oncs) > 0, "Empty class — check role strings above"

z_map = lambda ev: 'Z3' if ev <= -1.471 else ('Z1' if ev >= 0.382 else 'Z2')

# EXP1: Mann-Whitney continuous Ev
tsg_ev = [g['mean_ev'] for g in tsgs]
onc_ev = [g['mean_ev'] for g in oncs]
u1, p1 = mannwhitneyu(tsg_ev, onc_ev, alternative='less')
print(f"\nEXP1 Mann-Whitney TSG_ev < ONC_ev: U={u1:.0f}, p={p1:.3e}")
print(f"  TSG: mean={np.mean(tsg_ev):.3f} median={np.median(tsg_ev):.3f} std={np.std(tsg_ev):.3f}")
print(f"  ONC: mean={np.mean(onc_ev):.3f} median={np.median(onc_ev):.3f} std={np.std(onc_ev):.3f}")

# EXP2: Fisher Z3 vs Z1
z3_tsg = sum(1 for g in tsgs if g['mean_ev'] <= -1.471)
z1_tsg = sum(1 for g in tsgs if g['mean_ev'] >= 0.382)
z3_onc = sum(1 for g in oncs if g['mean_ev'] <= -1.471)
z1_onc = sum(1 for g in oncs if g['mean_ev'] >= 0.382)
ct2 = [[z3_tsg, z3_onc], [z1_tsg, z1_onc]]
or2, p2 = fisher_exact(ct2, alternative='greater')
print(f"\nEXP2 Fisher 2x2 [Z3_TSG={z3_tsg} Z3_ONC={z3_onc}] / [Z1_TSG={z1_tsg} Z1_ONC={z1_onc}]")
print(f"  OR={or2:.2f}, p={p2:.3e}")

# EXP3: Spearman Ev vs %Z3 within each class
r_tsg, p_tsg = spearmanr([g['mean_ev'] for g in tsgs], [g['pct_z3'] for g in tsgs])
r_onc, p_onc = spearmanr([g['mean_ev'] for g in oncs], [g['pct_z3'] for g in oncs])
print(f"\nEXP3 Spearman(mean_ev, pct_z3): TSG r={r_tsg:.3f} p={p_tsg:.3e} | ONC r={r_onc:.3f} p={p_onc:.3e}")

# EXP4: Chi2 full Z1/Z2/Z3 distribution
tsg_zones = [z_map(g['mean_ev']) for g in tsgs]
onc_zones = [z_map(g['mean_ev']) for g in oncs]
print(f"\nEXP4 Zone distribution:")
for z in ['Z1','Z2','Z3']:
    print(f"  {z}: TSG={tsg_zones.count(z)}/{len(tsgs)} ({tsg_zones.count(z)/len(tsgs):.1%})  "
          f"ONC={onc_zones.count(z)}/{len(oncs)} ({onc_zones.count(z)/len(oncs):.1%})")
ct4 = [[tsg_zones.count(z) for z in ['Z1','Z2','Z3']],
       [onc_zones.count(z) for z in ['Z1','Z2','Z3']]]
chi2, p4, dof, _ = chi2_contingency(ct4)
print(f"  Chi2={chi2:.2f}, p={p4:.3e}, dof={dof}")

# EXP5: n_windows-weighted Ev (heavier gene bodies more reliable)
tsg_wev = np.average([g['mean_ev'] for g in tsgs], weights=[g['n_windows'] for g in tsgs])
onc_wev = np.average([g['mean_ev'] for g in oncs], weights=[g['n_windows'] for g in oncs])
print(f"\nEXP5 n_windows-weighted mean_ev: TSG={tsg_wev:.3f}  ONC={onc_wev:.3f}  delta={tsg_wev-onc_wev:.3f}")

# EXP6: Extreme gene spotlight
print(f"\nEXP6 Most Z3-embedded TSGs:")
for g in sorted(tsgs, key=lambda x: x['mean_ev'])[:5]:
    print(f"  {g['gene']:12s} mean_ev={g['mean_ev']:.3f}  pct_z3={g['pct_z3']:.0%}  chr{g['chrom']}")
print(f"Most Z1-embedded Oncogenes:")
for g in sorted(oncs, key=lambda x: x['mean_ev'], reverse=True)[:5]:
    print(f"  {g['gene']:12s} mean_ev={g['mean_ev']:.3f}  pct_z3={g['pct_z3']:.0%}  chr{g['chrom']}")

results = {
    "exp1_mannwhitney": {"U": float(u1), "p": float(p1),
                          "tsg_mean": float(np.mean(tsg_ev)), "onc_mean": float(np.mean(onc_ev))},
    "exp2_fisher":      {"OR": float(or2), "p": float(p2), "contingency": ct2},
    "exp3_spearman":    {"tsg_r": float(r_tsg), "tsg_p": float(p_tsg),
                          "onc_r": float(r_onc), "onc_p": float(p_onc)},
    "exp4_chi2":        {"chi2": float(chi2), "p": float(p4),
                          "tsg_zones": {z: tsg_zones.count(z) for z in ['Z1','Z2','Z3']},
                          "onc_zones": {z: onc_zones.count(z) for z in ['Z1','Z2','Z3']}},
    "exp5_weighted_ev": {"tsg": float(tsg_wev), "onc": float(onc_wev)},
}
with open(OUT, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved: {OUT}")
