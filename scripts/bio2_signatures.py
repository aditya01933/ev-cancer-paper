#!/usr/bin/env python3
"""Mutation signature decomposition for LUAD, SKCM, PRAD.
Tests whether OR variation encodes mutational biology.
Expected: LUAD C→A Z1-enriched (tobacco), SKCM C→T Z3-enriched (UV),
PRAD identifies driver of highest OR=1.91."""
import urllib.request, gzip, json, os
import numpy as np
from scipy.stats import fisher_exact

P_CHECK = np.random.default_rng(42).standard_normal((256,500))/np.sqrt(256)
assert abs(P_CHECK[0,0]-0.01904482)<1e-5, "WRONG RNG"

CANCER_DIR=os.path.expanduser('~/ai-project/AncientKeyGen1/imp-research/tcga_cancer_experiment')
Z1_T, Z3_T = 0.382, -1.471
WS = 5000

# Load genome-wide zone lookup (all chromosomes)
print("Loading genome-wide zones...")
zone_lookup = {}
window_counts = {'z1':0,'z3':0}
for chrom in [str(i) for i in range(1,23)]+['X']:
    zf = os.path.join(CANCER_DIR,f'human_chr{chrom}_ev_zones.json')
    if not os.path.exists(zf): continue
    data = json.load(open(zf))
    zone_lookup[chrom] = {d['start']:d['ev_resid'] for d in data}
    window_counts['z1'] += sum(1 for d in data if d['ev_resid']>=Z1_T)
    window_counts['z3'] += sum(1 for d in data if d['ev_resid']<=Z3_T)
print(f"Loaded {sum(len(v) for v in zone_lookup.values()):,} windows")
print(f"Z1={window_counts['z1']:,}, Z3={window_counts['z3']:,}")

# Sanity
er_gene=zone_lookup.get('1',{}).get((5_000_000//WS)*WS)
er_peri=zone_lookup.get('1',{}).get((121_000_000//WS)*WS)
print(f"\nSanity: gene_rich={er_gene:+.3f}, pericentro={er_peri:+.3f}")

def get_zone(chrom, pos):
    w=(int(pos)//WS)*WS
    er=zone_lookup.get(str(chrom).replace('chr',''),{}).get(w)
    if er is None: return None
    return 'Z1' if er>=Z1_T else ('Z3' if er<=Z3_T else 'Z2')

def classify_mut(ref, alt):
    """Classify mutation to canonical pyrimidine-based SBS type"""
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    if ref in ('C','T'):
        return f"{ref}>{alt}"
    else:
        return f"{comp.get(ref,'?')}>{comp.get(alt,'?')}"

# === Process three cancer types ===
PROJECTS = {
    'TCGA_LUAD': {'label':'Lung adenocarcinoma', 'expected_high':'C>A (tobacco)'},
    'TCGA_SKCM': {'label':'Melanoma',            'expected_high':'C>T (UV)'},
    'TCGA_PRAD': {'label':'Prostate',             'expected_high':'unknown (OR=1.91)'},
}

all_results = {}

for project, meta in PROJECTS.items():
    print(f"\n{'='*55}")
    print(f"{project} — {meta['label']}")
    print(f"Expected: {meta['expected_high']}")
    print(f"{'='*55}")

    maf_cache = os.path.join(CANCER_DIR, f"maf_ids_{project}.json")
    if not os.path.exists(maf_cache):
        print(f"  No MAF cache found for {project} — skipping")
        continue
    maf_ids = json.load(open(maf_cache))
    print(f"  MAF files: {len(maf_ids)}")

    # Count mutations by zone AND mutation type
    mut_zone = {}  # (mut_type, zone) -> count
    zone_total = {'Z1':0,'Z2':0,'Z3':0}
    n_done = 0

    for entry in maf_ids[:50]:  # 50 samples per cancer type
        fid = entry['file_id'] if isinstance(entry,dict) else entry
        try:
            with urllib.request.urlopen(
                urllib.request.Request(
                    f"https://api.gdc.cancer.gov/data/{fid}",
                    headers={"User-Agent":"Mozilla/5.0"}),timeout=30) as r:
                raw = r.read()
            try: lines=gzip.decompress(raw).decode('utf-8','replace').split('\n')
            except: lines=raw.decode('utf-8','replace').split('\n')

            hdr=None; ds=0
            for i,l in enumerate(lines):
                if l.startswith('#'): continue
                hdr=l.split('\t'); ds=i+1; break
            if not hdr: continue

            cc=next((i for i,h in enumerate(hdr) if h=='Chromosome'),None)
            pc=next((i for i,h in enumerate(hdr) if h=='Start_Position'),None)
            rc=next((i for i,h in enumerate(hdr) if h=='Reference_Allele'),None)
            ac=next((i for i,h in enumerate(hdr) if h=='Tumor_Seq_Allele2'),None)
            vc=next((i for i,h in enumerate(hdr) if h=='Variant_Classification'),None)
            if None in (cc,pc,rc,ac): continue

            for line in lines[ds:]:
                if not line: continue
                cols=line.split('\t')
                if len(cols)<=max(c for c in [cc,pc,rc,ac] if c): continue
                try:
                    # SNPs only
                    ref=cols[rc]; alt=cols[ac]
                    if len(ref)!=1 or len(alt)!=1 or ref==alt: continue
                    if ref not in 'ACGT' or alt not in 'ACGT': continue
                    chrom=cols[cc]; pos=int(cols[pc])
                    z=get_zone(chrom,pos)
                    if not z: continue
                    mut=classify_mut(ref,alt)
                    mut_zone[(mut,z)]=mut_zone.get((mut,z),0)+1
                    zone_total[z]+=1
                except: continue
            n_done+=1
        except: continue

    print(f"  Processed: {n_done} MAFs")
    print(f"  Mutations: Z1={zone_total['Z1']:,}, Z3={zone_total['Z3']:,}")

    # Compute per-mutation-type Z3 vs Z1 enrichment
    mut_types = sorted(set(k[0] for k in mut_zone.keys()))
    print(f"\n  {'Mutation':<8} {'Z1_n':>7} {'Z3_n':>7} {'Z1%':>7} {'Z3%':>7} {'ratio':>7} {'p':>10} sig")
    print(f"  {'-'*65}")

    sig_results = []
    for mut in mut_types:
        z1n = mut_zone.get((mut,'Z1'),0)
        z3n = mut_zone.get((mut,'Z3'),0)
        if z1n+z3n < 20: continue
        z1pct = z1n/max(zone_total['Z1'],1)*100
        z3pct = z3n/max(zone_total['Z3'],1)*100
        ratio = (z3pct/z1pct) if z1pct>0 else 0
        # Fisher: is this mut type enriched in Z3 vs Z1?
        other_z1 = zone_total['Z1']-z1n
        other_z3 = zone_total['Z3']-z3n
        _,p=fisher_exact([[z3n,other_z3],[z1n,other_z1]])
        sig='🔴 Z3' if p<0.05 and ratio>1 else ('🔵 Z1' if p<0.05 and ratio<1 else '  ns')
        print(f"  {mut:<8} {z1n:>7,} {z3n:>7,} {z1pct:>6.1f}% {z3pct:>6.1f}% {ratio:>7.3f} {p:>10.2e} {sig}")
        sig_results.append({'mut':mut,'z1n':z1n,'z3n':z3n,
                            'z1pct':z1pct,'z3pct':z3pct,'ratio':ratio,'p':p})

    all_results[project] = {
        'label': meta['label'],
        'n_mafs': n_done,
        'zone_total': zone_total,
        'mutations': sig_results
    }

# === Summary across cancers ===
print(f"\n{'='*55}")
print("CROSS-CANCER SIGNATURE SUMMARY")
print(f"{'='*55}")
print("Expected biology:")
print("  LUAD: C>A enriched Z1 (tobacco=euchromatic)")
print("  SKCM: C>T enriched Z3 (UV at nucleosomes)")
print("  PRAD: identify driver of OR=1.91")
print()
for project, res in all_results.items():
    z3_enriched = [r['mut'] for r in res['mutations'] if r['p']<0.05 and r['ratio']>1]
    z1_enriched = [r['mut'] for r in res['mutations'] if r['p']<0.05 and r['ratio']<1]
    print(f"  {project}: Z3-enriched={z3_enriched}, Z1-enriched={z1_enriched}")

out_path = os.path.join(CANCER_DIR,'pan_cancer_signatures.json')
with open(out_path,'w') as f:
    json.dump(all_results,f,indent=2)
print(f"\n✅ Saved: pan_cancer_signatures.json")
