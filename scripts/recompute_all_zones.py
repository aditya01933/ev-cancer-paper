# scripts/recompute_all_zones.py
"""
Full recompute from raw data. Validates all paper claims.
Inputs:  raw_data/fasta/human_chrN.fa.gz
         raw_data/maf/mc3.v0.2.8.PUBLIC.maf.gz
Outputs: raw_zones/human_chrN_ev_zones_WITH_MUTATIONS.json
         data/recompute_summary.json
"""
import gzip, json, numpy as np, sys, os
from collections import defaultdict
from scipy.stats import skew as skewness, fisher_exact

sys.path.insert(0, 'scripts')
from ev_formula import P, NULL, KI, PRIMES, Z1_T, Z3_T

# ── 0. Formula fingerprint ──
assert abs(P[0,0] - 0.01904482) < 1e-6, "WRONG P MATRIX"
assert Z1_T == 0.382 and Z3_T == -1.471
assert len(PRIMES) == 95
print(f"[OK] Formula: P[0,0]={P[0,0]:.8f}  Z1={Z1_T}  Z3={Z3_T}")

TSS = {
    'A1':'BRCA','A2':'BRCA','A7':'BRCA','A8':'BRCA','AN':'BRCA','AO':'BRCA',
    'AQ':'BRCA','AR':'BRCA','B6':'BRCA','BH':'BRCA','C8':'BRCA','D8':'BRCA',
    'E2':'BRCA','E9':'BRCA','EW':'BRCA','GM':'BRCA','LL':'BRCA','OL':'BRCA',
    'PL':'BRCA','QS':'BRCA','S3':'BRCA','WT':'BRCA','XX':'BRCA',
    '05':'LUAD','38':'LUAD','49':'LUAD','4B':'LUAD','55':'LUAD','67':'LUAD',
    '6A':'LUAD','78':'LUAD','86':'LUAD','95':'LUAD','97':'LUAD','CG':'LUAD',
    'J2':'LUAD','MP':'LUAD','NJ':'LUAD','O1':'LUAD','QC':'LUAD','UB':'LUAD',
    '18':'LUSC','22':'LUSC','33':'LUSC','34':'LUSC','37':'LUSC','43':'LUSC',
    '56':'LUSC','60':'LUSC','63':'LUSC','66':'LUSC','77':'LUSC','85':'LUSC',
    '90':'LUSC','94':'LUSC','99':'LUSC','XF':'LUSC',
    'BF':'SKCM','D9':'SKCM','DA':'SKCM','EE':'SKCM','ER':'SKCM','FW':'SKCM',
    'GN':'SKCM','HR':'SKCM','IB':'SKCM','RP':'SKCM','WE':'SKCM','XK':'SKCM','YG':'SKCM',
    'CH':'PRAD','EJ':'PRAD','FC':'PRAD','G9':'PRAD','HC':'PRAD','KK':'PRAD',
    'M7':'PRAD','PD':'PRAD','QQ':'PRAD','V1':'PRAD','YL':'PRAD',
    '2F':'BLCA','4Z':'BLCA','BL':'BLCA','BT':'BLCA','DK':'BLCA','E7':'BLCA',
    'FD':'BLCA','GC':'BLCA','GU':'BLCA','HQ':'BLCA','K4':'BLCA','ZF':'BLCA',
    'A6':'COAD','AA':'COAD','AH':'COAD','AY':'COAD','D5':'COAD','DC':'COAD',
    'DM':'COAD','F4':'COAD','G4':'COAD','NH':'COAD','QG':'COAD','WS':'COAD',
    '02':'GBM','06':'GBM','08':'GBM','12':'GBM','14':'GBM','15':'GBM',
    '16':'GBM','19':'GBM','32':'GBM','41':'GBM','76':'GBM',
    '4P':'HNSC','BA':'HNSC','BB':'HNSC','BC':'HNSC','CN':'HNSC','CQ':'HNSC',
    'CV':'HNSC','D6':'HNSC','F7':'HNSC','HD':'HNSC','MT':'HNSC','UF':'HNSC',
    'A3':'KIRC','AK':'KIRC','B0':'KIRC','B2':'KIRC','B8':'KIRC','BP':'KIRC',
    'CJ':'KIRC','DZ':'KIRC','J9':'KIRC','MH':'KIRC','SX':'KIRC','Y8':'KIRC',
    'DD':'LIHC','ED':'LIHC','FV':'LIHC','G3':'LIHC','K7':'LIHC','MI':'LIHC',
    'RC':'LIHC','YA':'LIHC','ZS':'LIHC',
    '04':'OV','09':'OV','10':'OV','13':'OV','17':'OV','20':'OV','23':'OV',
    '24':'OV','25':'OV','26':'OV','27':'OV','29':'OV','30':'OV','31':'OV',
    '36':'OV','42':'OV','57':'OV','58':'OV','59':'OV','61':'OV','64':'OV','80':'OV',
    'BR':'STAD','D7':'STAD','HU':'STAD','IN':'STAD','RD':'STAD','W6':'STAD',
    'BJ':'THCA','BQ':'THCA','DY':'THCA','EL':'THCA','ET':'THCA','HF':'THCA',
    'J8':'THCA','KB':'THCA','QK':'THCA',
    'AP':'UCEC','AX':'UCEC','B5':'UCEC','BG':'UCEC','BS':'UCEC','D1':'UCEC',
    'D3':'UCEC','DF':'UCEC','DI':'UCEC','E6':'UCEC','EI':'UCEC','EO':'UCEC',
    'EY':'UCEC','FI':'UCEC','FS':'UCEC','N1':'UCEC',
}

CANCERS = ['BRCA','LUAD','LUSC','SKCM','PRAD','BLCA','COAD','GBM',
           'HNSC','KIRC','LIHC','OV','STAD','THCA','UCEC']
CHROMS  = [str(c) for c in range(1,23)] + ['X']
WS      = 5000
MAF     = 'raw_data/maf/mc3.v0.2.8.PUBLIC.maf.gz'
FASTA   = 'raw_data/fasta'
OUT     = 'raw_zones'

# ── 1. Parse MAF ──
print("\n[STEP 1] Parsing MAF...", flush=True)
mut = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
n_total, n_skip, n_mapped = 0, 0, 0

with gzip.open(MAF, 'rt') as f:
    header = None
    for line in f:
        if line.startswith('#'): continue
        row = line.strip().split('\t')
        if header is None:
            header = row
            ci = {k:i for i,k in enumerate(header)}
            assert 'Chromosome'           in ci
            assert 'Start_Position'       in ci
            assert 'Tumor_Sample_Barcode' in ci
            print(f"  MAF header OK: {len(header)} columns")
            continue
        try:
            chrom  = row[ci['Chromosome']].replace('chr','')
            pos    = int(row[ci['Start_Position']])
            tsb    = row[ci['Tumor_Sample_Barcode']]
            tss    = tsb.split('-')[1] if tsb.startswith('TCGA-') else ''
            cancer = TSS.get(tss, 'OTHER')
        except (IndexError, ValueError):
            n_skip += 1; continue
        if chrom not in set(CHROMS): continue
        mut[chrom][(pos // WS) * WS][cancer] += 1
        mut[chrom][(pos // WS) * WS]['ALL']   += 1
        n_total += 1
        if cancer != 'OTHER': n_mapped += 1

assert n_total  > 1_000_000, f"FAIL: {n_total} mutations"
assert n_mapped > 100_000,   f"FAIL: {n_mapped} mapped mutations"
assert len(mut) >= 22,       f"FAIL: {len(mut)} chroms"
print(f"  Total={n_total:,}  mapped={n_mapped:,} ({100*n_mapped/n_total:.1f}%)  skip={n_skip}")

cancer_totals = defaultdict(int)
for chrom in mut:
    for ws in mut[chrom]:
        for c,n in mut[chrom][ws].items():
            cancer_totals[c] += n
for c in CANCERS:
    assert cancer_totals[c] > 1000, f"FAIL: {c} only {cancer_totals[c]} mutations"
    print(f"  {c}: {cancer_totals[c]:,}")

# ── 2. Ev + GC helper ──
def ev_and_gc(seg):
    seg = seg.upper().replace('N','')
    if len(seg) < 1000: return None, None
    v = np.zeros(256)
    for i in range(len(seg)-3):
        k = seg[i:i+4]
        if k in KI: v[KI[k]] += 1
    if v.sum() == 0: return None, None
    f  = v / v.sum()
    gc = (seg.count('G') + seg.count('C')) / len(seg)
    e  = abs(skewness((f @ P)[PRIMES])) * 6.07 + 0.10
    er = (e - (NULL['slope'] * gc + NULL['intercept'])) / NULL['std']
    return round(float(er), 6), round(float(gc), 6)

# ── 3. Per-chromosome compute ──
print("\n[STEP 2] Computing zones...", flush=True)
summary = []

for chrom in CHROMS:
    fa = f'{FASTA}/human_chr{chrom}.fa.gz'
    assert os.path.exists(fa), f"MISSING: {fa}"

    with gzip.open(fa, 'rt') as f:
        seq = ''.join(l.strip() for l in f if not l.startswith('>'))
    assert len(seq) > 1_000_000, f"chr{chrom}: seq too short"

    windows, skipped = [], 0
    for start in range(0, len(seq) - WS, WS):
        er, gc = ev_and_gc(seq[start:start+WS])
        if er is None or not np.isfinite(er):
            skipped += 1; continue
        cm = dict(mut[chrom].get(start, {}))
        windows.append({
            'start':          start,
            'end':            start + WS,
            'ev_resid':       er,
            'gc':             gc,
            'mutation_count': int(cm.get('ALL', 0)),
            'cancer_counts':  {c: int(cm.get(c, 0)) for c in CANCERS}
        })

    ev_arr  = np.array([w['ev_resid']       for w in windows])
    gc_arr  = np.array([w['gc']             for w in windows])
    mut_arr = np.array([w['mutation_count'] for w in windows])
    n_z1 = int((ev_arr >= Z1_T).sum())
    n_z3 = int((ev_arr <= Z3_T).sum())

    assert len(windows) > 1000,      f"chr{chrom}: too few windows"
    assert n_z1 > 0,                 f"chr{chrom}: no Z1 windows"
    assert n_z3 > 0,                 f"chr{chrom}: no Z3 windows"
    assert ev_arr.std() > 0.1,       f"chr{chrom}: Ev variance too low"
    assert gc_arr.min() >= 0,        f"chr{chrom}: GC < 0"
    assert gc_arr.max() <= 1,        f"chr{chrom}: GC > 1"
    assert (mut_arr < 0).sum() == 0, f"chr{chrom}: negative mut counts"

    print(f"  chr{chrom:>2}: {len(windows):>6}w  skip={skipped:>4}  mut={mut_arr.sum():>8}  "
          f"Z1={n_z1:>5}  Z3={n_z3:>5}  gc=[{gc_arr.min():.3f},{gc_arr.max():.3f}]  "
          f"ev_mean={ev_arr.mean():+.3f}", flush=True)

    summary.append({'chrom': chrom, 'windows': len(windows),
                    'mutations': int(mut_arr.sum()), 'Z1': n_z1, 'Z3': n_z3})
    json.dump(windows, open(f'{OUT}/human_chr{chrom}_ev_zones_WITH_MUTATIONS.json', 'w'))

# ── 4. Genome-wide validation ──
print("\n[STEP 3] Genome-wide validation...", flush=True)

# Collect GC percentile thresholds first
all_gc_flat = []
for chrom in CHROMS:
    d = json.load(open(f'{OUT}/human_chr{chrom}_ev_zones_WITH_MUTATIONS.json'))
    all_gc_flat.extend(w['gc'] for w in d)
all_gc_arr = np.array(all_gc_flat)
p_z1_gc = np.percentile(all_gc_arr, 73.8)   # top 26.2% = Z1 proportion
p_z3_gc = np.percentile(all_gc_arr, 24.4)   # bottom 24.4% = Z3 proportion
print(f"  GC thresholds: Z1>={p_z1_gc:.4f}  Z3<={p_z3_gc:.4f}")

# Single pass: contingency tables for Ev-zones and GC-zones
# Use mutated-window (binary) for Fisher — avoids negative cell bug
ev_mw_z3, ev_uw_z3, ev_mw_z1, ev_uw_z1 = 0, 0, 0, 0
gc_mw_z3, gc_uw_z3, gc_mw_z1, gc_uw_z1 = 0, 0, 0, 0

for chrom in CHROMS:
    d = json.load(open(f'{OUT}/human_chr{chrom}_ev_zones_WITH_MUTATIONS.json'))
    for w in d:
        m  = w['mutation_count'] > 0   # binary: mutated or not
        ev = w['ev_resid']
        gc = w['gc']
        # Ev zones
        if ev >= Z1_T:
            if m: ev_mw_z1 += 1
            else: ev_uw_z1 += 1
        elif ev <= Z3_T:
            if m: ev_mw_z3 += 1
            else: ev_uw_z3 += 1
        # GC zones
        if gc >= p_z1_gc:
            if m: gc_mw_z1 += 1
            else: gc_uw_z1 += 1
        elif gc <= p_z3_gc:
            if m: gc_mw_z3 += 1
            else: gc_uw_z3 += 1

# Sanity
assert ev_mw_z1 > 0 and ev_mw_z3 > 0, "Zero mutated Ev windows"
assert ev_uw_z1 > 0 and ev_uw_z3 > 0, "Zero unmutated Ev windows"
assert gc_mw_z1 > 0 and gc_mw_z3 > 0, "Zero mutated GC windows"

ev_or, ev_p = fisher_exact([[ev_mw_z3, ev_uw_z3],[ev_mw_z1, ev_uw_z1]], alternative='greater')
gc_or, gc_p = fisher_exact([[gc_mw_z3, gc_uw_z3],[gc_mw_z1, gc_uw_z1]], alternative='greater')

total_w = sum(s['windows']   for s in summary)
total_m = sum(s['mutations'] for s in summary)

assert total_w > 500_000, f"FAIL: {total_w} windows"
assert total_m >  50_000, f"FAIL: {total_m} mutations"
assert 1.4 < ev_or < 2.2,  f"FAIL: Ev OR={ev_or:.3f} outside [1.4,2.2]"
assert ev_p < 1e-100,       f"FAIL: Ev p={ev_p:.2e}"

print(f"  Total windows : {total_w:,}")
print(f"  Total mutations: {total_m:,}")
print(f"\n  === KEY RESULTS ===")
print(f"  Ev-zones  OR={ev_or:.3f}  p={ev_p:.2e}  (paper claim: OR=1.682)")
print(f"  GC-zones  OR={gc_or:.3f}  p={gc_p:.2e}")
print(f"  Delta OR (Ev - GC) = {ev_or - gc_or:+.3f}")
print(f"  Ev OR replication: {'PASS' if abs(ev_or-1.682)<0.15 else 'DIVERGED — investigate'}")

result = {
    'total_windows': total_w, 'total_mutations': total_m,
    'ev_or': round(ev_or,4),  'ev_p': ev_p,
    'gc_or': round(gc_or,4),  'gc_p': gc_p,
    'delta_or': round(ev_or-gc_or, 4),
    'gc_z1_threshold': round(p_z1_gc,4),
    'gc_z3_threshold': round(p_z3_gc,4),
    'per_chrom': summary
}
json.dump(result, open('data/recompute_summary.json','w'), indent=2)
print(f"\n[DONE] → data/recompute_summary.json")
print("Verify: python3 scripts/sanity_zones.py")
