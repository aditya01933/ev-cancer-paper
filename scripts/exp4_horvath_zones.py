#!/usr/bin/env python3
"""
EXP 4: Horvath 353 Clock CpGs × Ev Zones
==========================================
Reads Horvath Additional File 3 CSV directly.
Maps CpG positions to your existing zone JSONs.
Tests: do clock CpGs cluster in specific Ev zones?
Splits by positive (193) vs negative (160) coefficient.

NOTE: Horvath positions are hg18 (Build 36).
      Zone JSONs are hg38. At 5kb resolution, most CpGs
      will map correctly despite ~few-hundred-bp offset.
      For rigorous analysis, liftOver to hg38 first.

Usage:
  python3 exp4_horvath_zones.py horvath_353_cpgs.csv /path/to/zone_jsons/
"""
import numpy as np, sys, json, csv, glob, os
from scipy.stats import chi2_contingency, mannwhitneyu, fisher_exact

Z1_T = +0.382
Z3_T = -1.471
WINDOW = 5000


def load_zones(zone_dir):
    """Load zone JSONs → {chr: list of window dicts}."""
    zones = {}
    for pattern in ['human_chr*_ev_zones.json', 'chr*_ev_zones.json']:
        files = sorted(glob.glob(os.path.join(zone_dir, pattern)))
        if files:
            break
    assert files, f"No zone JSONs in {zone_dir}"
    for f in files:
        bn = os.path.basename(f).replace('_ev_zones.json', '').replace('human_', '')
        chrom = bn if bn.startswith('chr') else 'chr' + bn
        with open(f) as fh:
            data = json.load(fh)
        if isinstance(data, list):
            zones[chrom] = data
        elif isinstance(data, dict):
            zones[chrom] = data.get('windows', list(data.values()))
        # Build position index for fast lookup
    total = sum(len(v) for v in zones.values())
    print(f"[ZONES] {len(zones)} chromosomes, {total:,} windows")
    return zones


def build_index(zones):
    """Build {chr: {window_start_bin: window_dict}} for O(1) lookup."""
    idx = {}
    for chrom, wins in zones.items():
        idx[chrom] = {}
        for w in wins:
            start = w.get('start', w.get('window_start', w.get('pos', None)))
            if start is None:
                continue
            binn = start // WINDOW
            idx[chrom][binn] = w
    return idx


def parse_horvath_csv(path):
    """Parse Horvath Additional File 3 → list of CpG dicts."""
    cpgs = []
    with open(path, encoding='utf-8', errors='replace') as f:
        # Skip header comment lines
        lines = f.readlines()

    # Find the actual header row (contains 'CpGmarker')
    header_idx = None
    for i, line in enumerate(lines):
        if 'CpGmarker' in line:
            header_idx = i
            break
    assert header_idx is not None, "Could not find 'CpGmarker' header row"

    import io
    # Auto-detect delimiter
    hdr_line = lines[header_idx]
    delim = '\t' if '\t' in hdr_line else ','
    reader = csv.DictReader(io.StringIO(''.join(lines[header_idx:])), delimiter=delim)
    for row in reader:
        cpg_id = row.get('CpGmarker', '').strip()
        if not cpg_id or cpg_id == '(Intercept)' or not cpg_id.startswith('cg'):
            continue
        try:
            chrom = 'chr' + row['Chr'].strip()
            pos = int(float(row['MapInfo'].strip()))
            coef = float(row['CoefficientTraining'].strip())
        except (ValueError, KeyError, TypeError):
            continue
        direction = 'positive' if coef > 0 else 'negative'
        # Prefer explicit column if available
        mar = row.get('Marginal Age Relationship', '').strip().lower()
        if mar in ('positive', 'negative'):
            direction = mar
        gene = row.get('Symbol', '').strip()
        cpgs.append({
            'cpg_id': cpg_id, 'chr': chrom, 'pos': pos,
            'coef': coef, 'direction': direction, 'gene': gene
        })
    print(f"[HORVATH] {len(cpgs)} CpGs parsed")
    n_pos = sum(1 for c in cpgs if c['direction'] == 'positive')
    n_neg = sum(1 for c in cpgs if c['direction'] == 'negative')
    print(f"  Positive (hypermethylated with age): {n_pos}")
    print(f"  Negative (hypomethylated with age):  {n_neg}")
    # ── Sanity: total should be ~353, split varies by definition ──
    assert n_pos + n_neg == len(cpgs), "pos+neg mismatch"
    assert len(cpgs) >= 340, f"Expected ~353 CpGs, got {len(cpgs)}"
    return cpgs


def map_cpgs_to_zones(cpgs, zone_idx):
    """Map each CpG to its zone window."""
    mapped = []
    missed = 0
    for cpg in cpgs:
        chrom = cpg['chr']
        if chrom not in zone_idx:
            missed += 1
            continue
        # Try exact bin and neighbors (handles hg18→hg38 offset)
        target_bin = cpg['pos'] // WINDOW
        w = None
        for offset in [0, -1, 1, -2, 2]:
            w = zone_idx[chrom].get(target_bin + offset)
            if w is not None:
                break
        if w is None:
            missed += 1
            continue
        ev_resid = w.get('ev_resid', w.get('evresid', None))
        if ev_resid is None:
            missed += 1
            continue
        zone = 'Z1' if ev_resid >= Z1_T else ('Z3' if ev_resid <= Z3_T else 'Z2')
        mapped.append({**cpg, 'ev_resid': float(ev_resid),
                       'gc': float(w.get('gc', 0)), 'zone': zone})
    print(f"[MAP] {len(mapped)} mapped, {missed} missed")
    assert len(mapped) >= 100, f"Too few mapped: {len(mapped)}. Check zone JSON format."
    return mapped


def genome_zone_distribution(zones):
    """Get genome-wide zone percentages for comparison."""
    counts = {'Z1': 0, 'Z2': 0, 'Z3': 0}
    total = 0
    for chrom, wins in zones.items():
        for w in wins:
            ev_resid = w.get('ev_resid', w.get('evresid', None))
            if ev_resid is None:
                continue
            total += 1
            if ev_resid >= Z1_T:
                counts['Z1'] += 1
            elif ev_resid <= Z3_T:
                counts['Z3'] += 1
            else:
                counts['Z2'] += 1
    pcts = {z: c / total * 100 for z, c in counts.items()}
    print(f"[GENOME] Z1={pcts['Z1']:.1f}% Z2={pcts['Z2']:.1f}% Z3={pcts['Z3']:.1f}% (n={total:,})")
    return counts, total


def test_enrichment(mapped, genome_counts, genome_total, label="ALL"):
    """Chi-squared + Fisher exact for zone enrichment."""
    cpg_counts = {'Z1': 0, 'Z2': 0, 'Z3': 0}
    for m in mapped:
        cpg_counts[m['zone']] += 1
    n = len(mapped)
    print(f"\n  [{label}] n={n}")
    print(f"  {'Zone':<6} {'Clock':>8} {'%':>8} {'Genome%':>10} {'Enrichment':>12}")
    for z in ['Z1', 'Z2', 'Z3']:
        pct_cpg = cpg_counts[z] / n * 100
        pct_genome = genome_counts[z] / genome_total * 100
        enrichment = pct_cpg / max(pct_genome, 0.01)
        print(f"  {z:<6} {cpg_counts[z]:>8} {pct_cpg:>7.1f}% {pct_genome:>9.1f}% {enrichment:>11.2f}x")

    # Chi-squared: clock vs genome distribution
    obs = np.array([cpg_counts['Z1'], cpg_counts['Z2'], cpg_counts['Z3']])
    exp_frac = np.array([genome_counts['Z1'], genome_counts['Z2'], genome_counts['Z3']]) / genome_total
    exp = exp_frac * n
    # ── Sanity: expected counts ──
    if (exp < 5).any():
        print(f"  [WARN] Expected counts < 5 in some cells: {exp}")
    chi2 = np.sum((obs - exp) ** 2 / exp)
    from scipy.stats import chi2 as chi2_dist
    p_chi2 = 1 - chi2_dist.cdf(chi2, df=2)
    print(f"  Chi-squared={chi2:.2f}, p={p_chi2:.2e}")

    # Fisher exact for Z1 specifically (expected enrichment for clock CpGs)
    z1_cpg = cpg_counts['Z1']
    non_z1_cpg = n - z1_cpg
    z1_genome = genome_counts['Z1']
    non_z1_genome = genome_total - z1_genome
    table_z1 = [[z1_cpg, non_z1_cpg], [z1_genome, non_z1_genome]]
    or_z1, p_z1 = fisher_exact(table_z1)
    print(f"  Z1 Fisher: OR={or_z1:.2f}, p={p_z1:.2e}")

    # Fisher exact for Z3
    z3_cpg = cpg_counts['Z3']
    non_z3_cpg = n - z3_cpg
    z3_genome = genome_counts['Z3']
    non_z3_genome = genome_total - z3_genome
    table_z3 = [[z3_cpg, non_z3_cpg], [z3_genome, non_z3_genome]]
    or_z3, p_z3 = fisher_exact(table_z3)
    print(f"  Z3 Fisher: OR={or_z3:.2f}, p={p_z3:.2e}")

    # ev_resid stats
    resids = [m['ev_resid'] for m in mapped]
    print(f"  Mean ev_resid={np.mean(resids):.4f} ± {np.std(resids):.4f}")

    return {
        'label': label, 'n': n, 'zone_counts': cpg_counts,
        'chi2': float(chi2), 'p_chi2': float(p_chi2),
        'z1_or': float(or_z1), 'z1_p': float(p_z1),
        'z3_or': float(or_z3), 'z3_p': float(p_z3),
        'mean_resid': float(np.mean(resids)),
        'std_resid': float(np.std(resids)),
    }


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 exp4_horvath_zones.py <horvath_353.csv> <zone_json_dir/>")
        return

    cpgs = parse_horvath_csv(sys.argv[1])
    zones = load_zones(sys.argv[2])
    zone_idx = build_index(zones)
    mapped = map_cpgs_to_zones(cpgs, zone_idx)
    genome_counts, genome_total = genome_zone_distribution(zones)

    print(f"\n{'='*60}")
    print(f"HORVATH 353 CLOCK CpGs × Ev ZONES")
    print(f"{'='*60}")

    # All 353
    r_all = test_enrichment(mapped, genome_counts, genome_total, "ALL 353")

    # Positive (hypermethylated with age) — expected: Z1 enriched (Polycomb/GC-rich)
    pos = [m for m in mapped if m['direction'] == 'positive']
    r_pos = test_enrichment(pos, genome_counts, genome_total, "POSITIVE (193, hypermeth)")

    # Negative (hypomethylated with age) — expected: Z2/Z3 shifted (shores/AT-rich)
    neg = [m for m in mapped if m['direction'] == 'negative']
    r_neg = test_enrichment(neg, genome_counts, genome_total, "NEGATIVE (160, hypometh)")

    # Compare positive vs negative ev_resid
    pos_resids = [m['ev_resid'] for m in pos]
    neg_resids = [m['ev_resid'] for m in neg]
    if len(pos_resids) >= 10 and len(neg_resids) >= 10:
        U, p_mw = mannwhitneyu(pos_resids, neg_resids, alternative='two-sided')
        print(f"\n  Positive vs Negative ev_resid: Mann-Whitney p={p_mw:.2e}")
        print(f"  Positive mean={np.mean(pos_resids):.4f}, Negative mean={np.mean(neg_resids):.4f}")
        diff = np.mean(pos_resids) - np.mean(neg_resids)
        print(f"  Difference={diff:.4f} ({'positive higher' if diff > 0 else 'negative higher'})")

    # ── Interpretation ──
    print(f"\n{'='*60}")
    print(f"INTERPRETATION")
    print(f"{'='*60}")
    if r_all['p_chi2'] > 0.05:
        print("  No significant zone enrichment for clock CpGs overall.")
    else:
        print(f"  Clock CpGs are NON-RANDOMLY distributed across Ev zones (p={r_all['p_chi2']:.2e})")

    if r_pos['z1_or'] > 1.5 and r_pos['z1_p'] < 0.05:
        print(f"  POSITIVE clock CpGs enriched in Z1 (OR={r_pos['z1_or']:.2f}, p={r_pos['z1_p']:.2e})")
        print(f"    → Consistent with Polycomb/PRC2 targets being GC-rich/euchromatic")
        print(f"    → Ev zones predict epigenetic aging clock sites from sequence alone")
    if r_neg['z3_or'] > 1.5 and r_neg['z3_p'] < 0.05:
        print(f"  NEGATIVE clock CpGs enriched in Z3 (OR={r_neg['z3_or']:.2f}, p={r_neg['z3_p']:.2e})")
        print(f"    → Consistent with CpG shores being AT-rich/heterochromatic")
    if r_pos['z1_or'] > 1.0 and r_neg['z1_or'] < r_pos['z1_or']:
        print(f"  Direction split confirmed: positive CpGs more Z1, negative CpGs less Z1")
        print(f"    → Aging clock has a COMPOSITIONAL basis detectable by Ev")

    # Top clock CpGs by |ev_resid|
    mapped_sorted = sorted(mapped, key=lambda m: abs(m['ev_resid']), reverse=True)
    print(f"\n  Top 10 clock CpGs by |ev_resid|:")
    print(f"  {'CpG':<14} {'Gene':<10} {'Dir':<6} {'Zone':<4} {'ev_resid':>10} {'coef':>10}")
    for m in mapped_sorted[:10]:
        print(f"  {m['cpg_id']:<14} {m['gene']:<10} {m['direction'][:3]:<6} "
              f"{m['zone']:<4} {m['ev_resid']:>10.3f} {m['coef']:>10.4f}")

    # Save
    out = {
        'all': r_all, 'positive': r_pos, 'negative': r_neg,
        'n_mapped': len(mapped), 'n_missed': len(cpgs) - len(mapped),
        'genome_zone_counts': genome_counts, 'genome_total': genome_total,
        'note': 'hg18 positions mapped to hg38 zones with ±2 bin tolerance',
    }
    with open('horvath_ev_zones_results.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: horvath_ev_zones_results.json")


if __name__ == '__main__':
    main()
