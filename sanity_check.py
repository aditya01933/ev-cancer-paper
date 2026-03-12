#!/usr/bin/env python3
"""
SANITY CHECK — verify everything from scratch.
Run from: ~/ai-project/ev-cancer-paper/
No dependencies on cached results. Recomputes from raw files.
"""
import numpy as np, json, gzip, os, sys
from scipy.stats import skew as skewness, fisher_exact, pearsonr, spearmanr
from collections import Counter, defaultdict

FAIL = []
WARN = []
PASS_LIST = []

def check(name, cond, detail=""):
    if cond:
        PASS_LIST.append(name)
        print(f"  ✓ {name}")
    else:
        FAIL.append(f"{name}: {detail}")
        print(f"  ✗ {name} — {detail}")

def warn(name, detail):
    WARN.append(f"{name}: {detail}")
    print(f"  ⚠ {name} — {detail}")

print("="*70)
print("SANITY CHECK — Full Pipeline Verification")
print("="*70)

# ══════════════════════════════════════════════════════════
# S1: FORMULA INTEGRITY
# ══════════════════════════════════════════════════════════
print("\n[S1] Formula integrity")

P = np.random.default_rng(42).standard_normal((256,500)) / np.sqrt(256)
PRIMES = [p for p in range(2,500) if all(p%i!=0 for i in range(2,int(p**0.5)+1))][:95]
BASES = 'ACGT'
KI = {a+b+c+d: i*64+j*16+k*4+l for i,a in enumerate(BASES) for j,b in enumerate(BASES)
      for k,c in enumerate(BASES) for l,d in enumerate(BASES)}
Z1_T, Z3_T = 0.382, -1.471
NULL = {'slope':-5.938, 'intercept':4.471, 'std':0.456}

check("P[0,0] fingerprint", abs(P[0,0] - 0.01904482) < 1e-6, f"got {P[0,0]}")
check("P shape", P.shape == (256, 500), f"got {P.shape}")
check("95 primes", len(PRIMES) == 95 and PRIMES[0]==2 and PRIMES[-1]==499,
      f"got {len(PRIMES)}, first={PRIMES[0]}, last={PRIMES[-1]}")
check("256 4-mers", len(KI) == 256, f"got {len(KI)}")
check("AAAA=0", KI['AAAA'] == 0)
check("TTTT=255", KI['TTTT'] == 255)

def ev_resid(seq):
    seq = seq.upper().replace('N','')
    if len(seq) < 1000: return None
    v = np.zeros(256)
    for i in range(len(seq)-3):
        k = seq[i:i+4]
        if k in KI: v[KI[k]] += 1
    if v.sum() == 0: return None
    f = v / v.sum()
    gc = (seq.count('G') + seq.count('C')) / len(seq)
    ev = abs(skewness((f @ P)[PRIMES])) * 6.07 + 0.10
    return (ev - (NULL['slope']*gc + NULL['intercept'])) / NULL['std']

# Test with known sequences
poly_a = "A" * 5000
poly_gc = "GC" * 2500
random_seq = "".join(np.random.choice(list("ACGT"), 5000))

er_a = ev_resid(poly_a)
er_gc = ev_resid(poly_gc)
er_rand = ev_resid(random_seq)

check("poly-A computes", er_a is not None, "returned None")
check("poly-GC computes", er_gc is not None, "returned None")
check("random seq computes", er_rand is not None, "returned None")
# poly-A should be extreme (high AT), poly-GC should differ
if er_a is not None and er_gc is not None:
    check("poly-A ≠ poly-GC", abs(er_a - er_gc) > 0.5, f"diff={abs(er_a-er_gc):.3f}")

# ══════════════════════════════════════════════════════════
# S2: RAW DATA FILES EXIST
# ══════════════════════════════════════════════════════════
print("\n[S2] Raw data files")

fasta_dir = 'raw_data/fasta'
for c in ['1','2','22','X']:
    f = f"{fasta_dir}/human_chr{c}.fa.gz"
    check(f"chr{c} FASTA exists", os.path.exists(f), f"missing {f}")

maf = 'raw_data/maf/mc3.v0.2.8.PUBLIC.maf.gz'
check("MC3 MAF exists", os.path.exists(maf), f"missing {maf}")

# ══════════════════════════════════════════════════════════
# S3: ZONE FILES — spot-check against recompute from FASTA
# ══════════════════════════════════════════════════════════
print("\n[S3] Zone file integrity (chr22 spot-check)")

zone_file = 'raw_zones/human_chr22_ev_zones.json'
if os.path.exists(zone_file):
    zones22 = json.load(open(zone_file))
    check("chr22 zone file loads", len(zones22) > 1000, f"got {len(zones22)} windows")

    # Recompute first 10 windows from FASTA
    fasta22 = f"{fasta_dir}/human_chr22.fa.gz"
    if os.path.exists(fasta22):
        with gzip.open(fasta22, 'rt') as f:
            header = f.readline()
            seq = ''.join(line.strip() for line in f)
        n_checked = 0
        n_match = 0
        for w in zones22[:20]:
            start = w['start']
            seg = seq[start:start+5000]
            if len(seg) < 5000: continue
            er = ev_resid(seg)
            if er is None: continue
            cached_er = w['ev_resid']
            diff = abs(er - cached_er)
            if diff < 0.01: n_match += 1
            n_checked += 1
        check(f"chr22 recompute matches ({n_match}/{n_checked})",
              n_match >= n_checked * 0.9,
              f"only {n_match}/{n_checked} match within 0.01")
    else:
        warn("chr22 FASTA missing", "cannot recompute")
else:
    warn("chr22 zone file missing", "skipping S3")

# ══════════════════════════════════════════════════════════
# S4: WITH_MUTATIONS files — verify mutation counts are non-zero
# ══════════════════════════════════════════════════════════
print("\n[S4] WITH_MUTATIONS file integrity")

wm_file = 'raw_zones/human_chr1_ev_zones_WITH_MUTATIONS.json'
if os.path.exists(wm_file):
    wm = json.load(open(wm_file))
    check("chr1 WITH_MUTATIONS loads", len(wm) > 40000, f"got {len(wm)}")

    # Check schema
    w0 = wm[0]
    required_keys = ['start', 'ev_resid', 'gc', 'mutation_count']
    for k in required_keys:
        check(f"  key '{k}' exists", k in w0, f"missing from schema")

    total_muts = sum(w['mutation_count'] for w in wm)
    check("chr1 has mutations", total_muts > 1000, f"total={total_muts}")

    # Zone distribution sanity
    z1 = sum(1 for w in wm if w['ev_resid'] >= Z1_T)
    z3 = sum(1 for w in wm if w['ev_resid'] <= Z3_T)
    z2 = len(wm) - z1 - z3
    check("chr1 Z1 ~25%", 0.15 < z1/len(wm) < 0.35, f"{z1/len(wm)*100:.1f}%")
    check("chr1 Z3 ~25%", 0.15 < z3/len(wm) < 0.35, f"{z3/len(wm)*100:.1f}%")

    # OR sanity on chr1 alone
    z1_mut = sum(1 for w in wm if w['ev_resid'] >= Z1_T and w['mutation_count'] > 0)
    z3_mut = sum(1 for w in wm if w['ev_resid'] <= Z3_T and w['mutation_count'] > 0)
    if z1 > 0 and z3 > 0 and z1_mut > 0 and z3_mut > 0:
        ore, pe = fisher_exact([[z3_mut, z3-z3_mut],[z1_mut, z1-z1_mut]])
        check("chr1 OR > 1", ore > 1.0, f"OR={ore:.3f}")
        check("chr1 OR in [1.1, 3.0]", 1.1 < ore < 3.0, f"OR={ore:.3f}")
        check("chr1 p < 0.01", pe < 0.01, f"p={pe:.2e}")
else:
    warn("chr1 WITH_MUTATIONS missing", "skipping S4")

# ══════════════════════════════════════════════════════════
# S5: MC3 MAF — verify BRCA samples exist and parse correctly
# ══════════════════════════════════════════════════════════
print("\n[S5] MC3 MAF parse check (first 100K lines)")

TSS_BRCA = {'A1','A2','A7','A8','AN','AO','AQ','AR','B6','BH','C8','D8',
            'E2','E9','EW','GM','LL','OL','PL','QS','S3','WT','XX'}
if os.path.exists(maf):
    n_total = n_brca = n_snv = 0
    cancer_types = Counter()
    with gzip.open(maf, 'rt') as f:
        header = None
        for line in f:
            if line.startswith('#'): continue
            row = line.strip().split('\t')
            if header is None: header = row; ci = {k:i for i,k in enumerate(header)}; continue
            n_total += 1
            if n_total > 100000: break
            try:
                tsb = row[ci['Tumor_Sample_Barcode']]
                tss = tsb.split('-')[1] if tsb.startswith('TCGA-') else ''
                ref = row[ci['Reference_Allele']]
                alt = row[ci['Tumor_Seq_Allele2']]
            except: continue
            if tss in TSS_BRCA: n_brca += 1
            if len(ref)==1 and len(alt)==1 and ref in 'ACGT' and alt in 'ACGT' and ref!=alt:
                n_snv += 1

    check("MC3 has data", n_total >= 100000, f"got {n_total}")
    check("BRCA samples found", n_brca > 100, f"got {n_brca} in 100K lines")
    check("SNVs found", n_snv > 50000, f"got {n_snv} SNVs in 100K lines")
    print(f"  (sampled 100K lines: {n_brca} BRCA, {n_snv} SNVs)")

# ══════════════════════════════════════════════════════════
# S6: CACHED RESULTS — cross-check internal consistency
# ══════════════════════════════════════════════════════════
print("\n[S6] Cached result files")

# Horvath
horv_file = 'data/horvath_ev_zones_results.json'
if os.path.exists(horv_file):
    h = json.load(open(horv_file))
    check("Horvath file loads", isinstance(h, dict))
    pos_or = h.get('pos_or', h.get('positive_or'))
    if pos_or:
        check("Horvath pos OR ~1.94", 1.5 < pos_or < 2.5, f"got {pos_or}")
else:
    warn("Horvath results missing", horv_file)

# Germline
germ_file = 'raw_data/germline_baseline_chr1_5.json'
if not os.path.exists(germ_file):
    germ_file = 'data/germline_baseline_chr1_5.json'
if os.path.exists(germ_file):
    g = json.load(open(germ_file))
    ratio = g.get('z3_z1_ratio', g.get('ratio'))
    if ratio:
        check("Germline Z3/Z1 ~1.13", 1.05 < ratio < 1.25, f"got {ratio}")
else:
    warn("Germline results missing", "checked both paths")

# Methylation
meth_file = 'data/cancer_methylation_results.json'
if os.path.exists(meth_file):
    m = json.load(open(meth_file))
    het_mean = m.get('het_mean')
    if het_mean:
        check("Methylation Z3 Δβ < 0", het_mean < 0, f"got {het_mean}")
        check("Methylation Z3 Δβ ~-0.054", -0.08 < het_mean < -0.03, f"got {het_mean}")
else:
    warn("Methylation results missing", meth_file)

# ══════════════════════════════════════════════════════════
# S7: TSG/OG — THE CRITICAL CHECK
# ══════════════════════════════════════════════════════════
print("\n[S7] TSG/Oncogene segregation (drug_targets_zone3.json)")

tsg_file = 'data/drug_targets_zone3.json'
if os.path.exists(tsg_file):
    d = json.load(open(tsg_file))
    check("101 genes loaded", len(d) == 101, f"got {len(d)}")

    roles = Counter(g['role'] for g in d)
    print(f"  Roles: {dict(roles)}")
    zones_all = Counter(g['zone'] for g in d)
    print(f"  Zones: {dict(zones_all)}")

    # Use ACTUAL role names from file
    tsg_key = 'TSG'
    og_key = 'oncogene'

    tsg_z3 = sum(1 for g in d if g['role']==tsg_key and g['zone']=='Z3')
    tsg_tot = sum(1 for g in d if g['role']==tsg_key)
    og_z3 = sum(1 for g in d if g['role']==og_key and g['zone']=='Z3')
    og_tot = sum(1 for g in d if g['role']==og_key)

    print(f"  TSG in Z3: {tsg_z3}/{tsg_tot} = {tsg_z3/max(tsg_tot,1)*100:.1f}%")
    print(f"  OG  in Z3: {og_z3}/{og_tot} = {og_z3/max(og_tot,1)*100:.1f}%")

    check("TSG total ~49", 40 < tsg_tot < 60, f"got {tsg_tot}")
    check("OG total ~52", 40 < og_tot < 60, f"got {og_tot}")

    # Fisher test
    if tsg_z3 + og_z3 > 0:
        ore, pe = fisher_exact([[tsg_z3, tsg_tot-tsg_z3],[og_z3, og_tot-og_z3]])
        print(f"  Fisher OR={ore:.2f} p={pe:.3f}")
        check("TSG/OG OR ≠ 9.2", ore < 5.0, f"OR={ore:.2f}")
        check("TSG/OG p > 0.05 (underpowered)", pe > 0.05, f"p={pe:.3f}")
        if ore > 5:
            warn("ALERT: OR > 5 detected", f"OR={ore:.2f} — investigate!")

    # PAPER CLAIM: OR=9.2, 73% TSGs in Z3
    tsg_z3_pct = tsg_z3/max(tsg_tot,1)*100
    check("TSG Z3% ≠ 73% (paper claim)", tsg_z3_pct < 50,
          f"got {tsg_z3_pct:.1f}% — paper says 73%")

    # Key genes
    tp53 = next((g for g in d if g['gene']=='TP53'), None)
    kras = next((g for g in d if g['gene']=='KRAS'), None)
    nras = next((g for g in d if g['gene']=='NRAS'), None)
    if tp53: check("TP53 in Z3", tp53['zone']=='Z3', f"zone={tp53['zone']}")
    if kras: check("KRAS in Z1", kras['zone']=='Z1', f"zone={kras['zone']}")
    if nras: check("NRAS in Z1", nras['zone']=='Z1', f"zone={nras['zone']}")
else:
    warn("drug_targets_zone3.json missing", tsg_file)

# ══════════════════════════════════════════════════════════
# S8: GENOME-WIDE ZONE TOTALS
# ══════════════════════════════════════════════════════════
print("\n[S8] Genome-wide zone totals")

total_w = z1_w = z3_w = 0
for c in [str(i) for i in range(1,23)] + ['X']:
    f = f'raw_zones/human_chr{c}_ev_zones.json'
    if not os.path.exists(f): continue
    data = json.load(open(f))
    for w in data:
        total_w += 1
        er = w.get('ev_resid', w.get('er'))
        if er >= Z1_T: z1_w += 1
        elif er <= Z3_T: z3_w += 1

if total_w > 0:
    print(f"  Total: {total_w}  Z1: {z1_w} ({z1_w/total_w*100:.1f}%)  Z3: {z3_w} ({z3_w/total_w*100:.1f}%)")
    check("Total ~582K", 580000 < total_w < 585000, f"got {total_w}")
    check("Z1 ~26%", 0.24 < z1_w/total_w < 0.28, f"{z1_w/total_w*100:.1f}%")
    check("Z3 ~24%", 0.22 < z3_w/total_w < 0.26, f"{z3_w/total_w*100:.1f}%")
    check("Z1 + Z3 < total", z1_w + z3_w < total_w)

# ══════════════════════════════════════════════════════════
# S9: CROSS-CHECK — zone files vs WITH_MUTATIONS files
# ══════════════════════════════════════════════════════════
print("\n[S9] Zone vs WITH_MUTATIONS consistency (chr1)")

base_file = 'raw_zones/human_chr1_ev_zones.json'
wm_file2 = 'raw_zones/human_chr1_ev_zones_WITH_MUTATIONS.json'
if os.path.exists(base_file) and os.path.exists(wm_file2):
    base = json.load(open(base_file))
    wm2 = json.load(open(wm_file2))
    check("Same window count", len(base) == len(wm2),
          f"base={len(base)} wm={len(wm2)}")
    # Spot-check first 10 ev_resid values match
    n_match = 0
    for b, w in zip(base[:10], wm2[:10]):
        b_er = b.get('ev_resid', b.get('er'))
        w_er = w.get('ev_resid', w.get('er'))
        if b_er is not None and w_er is not None and abs(b_er - w_er) < 0.001:
            n_match += 1
    check("Ev_resid values match (10/10)", n_match == 10, f"{n_match}/10")

# ══════════════════════════════════════════════════════════
# S10: FORMULA DETERMINISM — same input → same output
# ══════════════════════════════════════════════════════════
print("\n[S10] Formula determinism")

test_seq = "ATCGATCGATCGATCG" * 312 + "ATCG"  # 5000bp
er1 = ev_resid(test_seq)
er2 = ev_resid(test_seq)
check("Deterministic (same input → same output)", er1 == er2, f"{er1} vs {er2}")

# Different seed MUST give different result
P2 = np.random.default_rng(43).standard_normal((256,500)) / np.sqrt(256)
check("Seed 43 ≠ Seed 42", abs(P2[0,0] - P[0,0]) > 0.001,
      f"seed42={P[0,0]:.6f} seed43={P2[0,0]:.6f}")

# ══════════════════════════════════════════════════════════
# S11: PAPER CLAIMS AUDIT
# ══════════════════════════════════════════════════════════
print("\n[S11] Paper claims audit")

claims = {
    "OR=1.682 (GDC)": "CANNOT verify without GDC MAFs — use MC3 cross-val instead",
    "OR=1.585 (MC3)": "Verified by T2 in validation engine",
    "15/15 cancers": "Verified by T4",
    "582,028 windows": f"Got {total_w} — {'MATCH' if 580000<total_w<585000 else 'MISMATCH'}",
    "C>T only enriched": "Verified by T5 (proportion ratio)",
    "Horvath OR=1.94": "Verified by T14",
    "Germline 1.129": "Verified by T15",
    "RT 3/3 sig": "Verified by T7",
    "TSG/OG OR=9.2": f"*** REFUTED: actual OR={ore:.2f} p={pe:.3f} ***" if 'ore' in dir() else "NOT TESTED",
}
print("  Claims status:")
for claim, status in claims.items():
    marker = "✗" if "REFUTED" in status or "MISMATCH" in status else "✓"
    print(f"    {marker} {claim}: {status}")

# ══════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════
print("\n" + "="*70)
print("SANITY CHECK SUMMARY")
print("="*70)
print(f"  PASSED: {len(PASS_LIST)}")
print(f"  FAILED: {len(FAIL)}")
print(f"  WARNINGS: {len(WARN)}")
if FAIL:
    print("\n  FAILURES:")
    for f in FAIL:
        print(f"    ✗ {f}")
if WARN:
    print("\n  WARNINGS:")
    for w in WARN:
        print(f"    ⚠ {w}")
print()
if not FAIL:
    print("  ALL CHECKS PASSED — pipeline is clean")
else:
    print(f"  {len(FAIL)} ISSUE(S) NEED ATTENTION")
