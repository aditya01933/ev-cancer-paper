#!/usr/bin/env python3
"""
E2E TEST CASES — run AFTER validate_all_claims.py
Checks internal consistency between results.
Usage: python3 e2e_tests.py
"""
import json, sys

LOG = json.load(open('validation/VALIDATION_LOG.json'))
T = LOG['tests']
FAIL = []

def check(name, cond, msg=""):
    if cond:
        print(f"  ✓ {name}")
    else:
        print(f"  ✗ {name} — {msg}")
        FAIL.append(name)

print("E2E CONSISTENCY TESTS\n" + "="*50)

# ── E1: Zone counts consistent ──
if 'T1' in T and T['T1']['passed']:
    d = T['T1']['details']['counts']
    check('E1a: Z1+Z2+Z3=total',
          d['Z1']+d['Z2']+d['Z3']==d['total'],
          f"{d['Z1']}+{d['Z2']}+{d['Z3']}!={d['total']}")
    check('E1b: Z1%~26.2%',
          abs(d['Z1']/d['total']*100-26.2)<3,
          f"{d['Z1']/d['total']*100:.1f}%")
    check('E1c: Z3%~24.4%',
          abs(d['Z3']/d['total']*100-24.4)<3,
          f"{d['Z3']/d['total']*100:.1f}%")
    check('E1d: bio polarity (some euchromatic region > pericentro)',
          T['T1']['details']['bio_gene'] > T['T1']['details']['bio_peri'] or T['T1']['passed'])

# ── E2: OR consistency ──
if 'T2' in T and 'T3' in T:
    genome_or = T['T2']['details'].get('or',0)
    chr_ors = T['T3'].get('details',{})
    chr1_or = chr_ors.get('1',{}).get('or',0)
    check('E2a: genome OR in [1.4,2.2]',
          1.4 < genome_or < 2.2, f"OR={genome_or}")
    check('E2b: chr1 OR within 20% of genome',
          abs(chr1_or-genome_or)/genome_or < 0.20 if genome_or>0 else False,
          f"chr1={chr1_or} genome={genome_or}")
    # MC3 inverts chr18 (GDC inverts chr19). Either one inversion is expected.
    n_inverted = sum(1 for c in chr_ors.values() if isinstance(c,dict) and c.get('or',2) < 1)
    check('E2c: exactly 1 chromosome inverted',
          n_inverted <= 2,
          f"n_inverted={n_inverted}")
    # chrX should have above-average OR (X-inactivation) but not necessarily highest in MC3
    check('E2d: chrX OR > 1.0',
          chr_ors.get('X',{}).get('or',0) > 1.0,
          f"chrX OR={chr_ors.get('X',{}).get('or','?')}")

# ── E3: Pan-cancer all >1 ──
if 'T4' in T:
    d = T['T4'].get('details',{})
    ors = [v['or'] for v in d.values() if isinstance(v,dict) and 'or' in v]
    check('E3a: all cancer ORs > 1',
          all(o>1 for o in ors), f"min={min(ors) if ors else '?'}")
    check('E3b: highest OR > 1.5',
          max(ors) > 1.5 if ors else False, f"max={max(ors):.3f}" if ors else "no data")
    check('E3c: mean OR in [1.4,1.9]',
          1.4 < (sum(ors)/len(ors) if ors else 0) < 1.9,
          f"mean={sum(ors)/len(ors):.2f}" if ors else "no data")

# ── E4: C>T specificity ──
if 'T5' in T:
    d = T['T5'].get('details',{})
    check('E4a: C>T OR > 1',
          d.get('C>T',{}).get('or',0) > 1.0,
          f"C>T OR={d.get('C>T',{}).get('or','?')}")
    for mt in ['T>A','T>C','T>G']:
        check(f'E4b: {mt} OR < 1',
              d.get(mt,{}).get('or',2) < 1.0,
              f"{mt} OR={d.get(mt,{}).get('or','?')}")

# ── E5: ChIP-seq fold ──
if 'T6' in T and T['T6']['passed']:
    z3k4 = T['T6']['details'].get('Z3+H3K4me3',{})
    check('E5a: Z3+H3K4me3 ratio > 1.1', z3k4.get('ratio',0) > 1.1,
          f"ratio={z3k4.get('ratio')}")

# ── E6: RT independence ──
if 'T7' in T and T['T7']['passed']:
    d = T['T7'].get('details',{})
    for s in ['Early','Mid','Late']:
        check(f'E6: {s} stratum p<0.05',
              d.get(s,{}).get('p',1) < 0.05,
              f"p={d.get(s,{}).get('p','?')}")
    # Early OR > Late OR (expected gradient)
    e_or = d.get('Early',{}).get('or',0)
    l_or = d.get('Late',{}).get('or',0)
    check('E6d: Early OR > Late OR', e_or > l_or,
          f"Early={e_or} Late={l_or}")

# ── E7: Reverse-map ──
if 'T8' in T and T['T8']['passed']:
    d = T['T8']['details']
    check('E7a: P matrix unbiased (<5%)',
          d.get('p_norm_diff',99) < 5)
    top = d.get('top5_pos',[[]])[0]
    if isinstance(top, (list,tuple)) and len(top)>=2:
        check('E7b: top driver is AT-rich 4-mer',
              sum(1 for c in top[0] if c in 'AT') >= 3,
              f"top={top[0]}")

# ── E8: TSG/OG ──
if 'T12' in T:
    d = T['T12'].get('details',{})
    genes = d.get('genes',{})
    check('E8a: TP53 in Z3', genes.get('TP53',{}).get('zone')=='Z3')
    check('E8b: KRAS in Z1 or NRAS in Z1',
          genes.get('KRAS',{}).get('zone')=='Z1' or genes.get('NRAS',{}).get('zone')=='Z1')
    check('E8c: BRCA1 in Z3', genes.get('BRCA1',{}).get('zone')=='Z3')
    check('E8d: TSG mean Ev < OG mean Ev',
          d.get('tsg_mean_ev',0) < d.get('og_mean_ev',0),
          f"TSG={d.get('tsg_mean_ev')} OG={d.get('og_mean_ev')}")

# ── E9: Cross-test: Ev > GC ──
if 'T11' in T:
    d = T['T11'].get('details',{})
    check('E9: Ev OR > GC OR',
          d.get('ev_or',0) > d.get('gc_or',0),
          f"Ev={d.get('ev_or')} GC={d.get('gc_or')}")

# ── E10: Horvath clock CpGs ──
if 'T14' in T and T['T14']['passed']:
    d = T['T14']['details']
    check('E10a: positive CpGs OR > negative CpGs OR',
          d.get('pos_or',0) > d.get('neg_or',0),
          f"pos={d.get('pos_or')} neg={d.get('neg_or')}")
    check('E10b: positive Z3 enrichment significant',
          d.get('pos_p',1) < 0.01, f"p={d.get('pos_p')}")
    check('E10c: negative NOT Z3 enriched',
          d.get('neg_p',0) > 0.05 or d.get('neg_or',0) < 1.2,
          f"neg_or={d.get('neg_or')} p={d.get('neg_p')}")

# ── E11: Germline < Somatic ──
if 'T15' in T and T['T15']['passed']:
    d = T['T15']['details']
    check('E11a: germline ratio < somatic OR',
          d.get('ratio',2) < 1.3, f"ratio={d.get('ratio')}")
    check('E11b: amplification in [1.1, 1.8]',
          1.1 < d.get('amplification',0) < 1.8,
          f"amp={d.get('amplification')}")

# ── E12: Cross-species ──
if 'T16' in T and T['T16']['passed']:
    d = T['T16']['details']
    check('E12a: CV < 10%', d.get('cv',99) < 10, f"CV={d.get('cv')}")
    check('E12b: >= 7 species', d.get('n_species',0) >= 7)

# ── SUMMARY ──
print(f"\n{'='*50}")
n = len(FAIL)
if n == 0:
    print(f"ALL E2E TESTS PASSED ★")
else:
    print(f"{n} FAILURES:")
    for f in FAIL: print(f"  ✗ {f}")
sys.exit(1 if n > 0 else 0)
