# mut_signature_zones.py
# Which mutation types are enriched in Zone3 vs Zone1?
# Reuses cached MAF IDs and zone files

import urllib.request, json, gzip, os
import numpy as np
from scipy import stats
from collections import defaultdict

CANCER_DIR = os.path.expanduser("~/ai-project/AncientKeyGen1/imp-research/tcga_cancer_experiment")
os.chdir(CANCER_DIR)
Z3_T, Z1_T, WINDOW = -1.471, 0.382, 5000

# Load zones
print("Loading zones...")
zone_lookup = {}
for c in [str(i) for i in range(1,23)] + ['X']:
    zf = f"human_chr{c}_ev_zones.json"
    if not os.path.exists(zf): continue
    with open(zf) as f: zones = json.load(f)
    zone_lookup[c] = {w['start']: w['ev_resid'] for w in zones}

def get_zone(chrom, pos):
    c = chrom.replace('chr','')
    e = zone_lookup.get(c,{}).get((int(pos)//WINDOW)*WINDOW)
    if e is None: return None
    return 'Z3' if e<=Z3_T else ('Z1' if e>=Z1_T else 'Z2')

# Load MAF IDs — use BRCA (largest, already cached)
with open('maf_ids_TCGA_BRCA.json') as f: hits = json.load(f)
print(f"Using {len(hits)} BRCA MAFs")

# Count mutation types per zone
zone_muts = {'Z1': defaultdict(int), 'Z2': defaultdict(int), 'Z3': defaultdict(int)}
n_done = 0

for hit in hits[:200]:
    try:
        url = f"https://api.gdc.cancer.gov/data/{hit['file_id']}"
        req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as r: raw = r.read()
        try:    content = gzip.decompress(raw).decode('utf-8', errors='replace')
        except: content = raw.decode('utf-8', errors='replace')
        lines = content.split('\n')
        header, ds = None, 0
        for i, line in enumerate(lines):
            if line.startswith('#'): continue
            header = line.split('\t'); ds = i+1; break
        if not header: continue

        chr_c = next((i for i,h in enumerate(header) if h=='Chromosome'), None)
        pos_c = next((i for i,h in enumerate(header) if h=='Start_Position'), None)
        ref_c = next((i for i,h in enumerate(header) if h=='Reference_Allele'), None)
        alt_c = next((i for i,h in enumerate(header) if h=='Tumor_Seq_Allele2'), None)
        if any(x is None for x in [chr_c,pos_c,ref_c,alt_c]): continue

        for line in lines[ds:]:
            if not line: continue
            p = line.split('\t')
            if len(p) <= max(chr_c,pos_c,ref_c,alt_c): continue
            z = get_zone(p[chr_c], p[pos_c])
            if not z: continue
            ref, alt = p[ref_c], p[alt_c]
            if len(ref)==1 and len(alt)==1 and ref!='-' and alt!='-':
                zone_muts[z][f"{ref}>{alt}"] += 1
        n_done += 1
        if n_done % 50 == 0: print(f"  {n_done}/200")
    except: continue

print(f"\nProcessed {n_done} samples")

# Normalize and compare Z3 vs Z1
z1_total = sum(zone_muts['Z1'].values())
z3_total = sum(zone_muts['Z3'].values())
all_types = sorted(set(list(zone_muts['Z1'].keys()) + list(zone_muts['Z3'].keys())))

print(f"\n{'Mutation':<8} {'Z1%':>7} {'Z3%':>7} {'Z3/Z1':>7} {'enriched'}")
print("-"*40)

results = []
for mt in all_types:
    z1n = zone_muts['Z1'].get(mt,0)
    z3n = zone_muts['Z3'].get(mt,0)
    if z1n < 10: continue
    z1p = z1n/z1_total
    z3p = z3n/z3_total
    ratio = z3p/z1p if z1p>0 else 0
    ct = np.array([[z3n, z3_total-z3n],[z1n, z1_total-z1n]])
    _, p = stats.fisher_exact(ct, alternative='greater' if ratio>1 else 'less')
    results.append((mt, z1p, z3p, ratio, p))

results.sort(key=lambda x: -x[3])
for mt, z1p, z3p, ratio, p in results:
    flag = '🔴 Z3-enriched' if ratio>1.1 and p<0.05 else ('🔵 Z1-enriched' if ratio<0.9 and p<0.05 else '')
    print(f"{mt:<8} {z1p*100:>6.2f}% {z3p*100:>6.2f}% {ratio:>7.2f}x  {flag}")

json.dump([{'type':r[0],'z1_pct':r[1],'z3_pct':r[2],'ratio':r[3],'p':r[4]}
           for r in results], open('mutation_signature_zones.json','w'), indent=2)
print("\nSaved → mutation_signature_zones.json")
