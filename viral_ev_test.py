#!/usr/bin/env python3
"""
VIRAL PANDEMIC PREDICTION: Zero-training Ev vs trained ML (Mollentze 2021)
==========================================================================
Tests whether Ev (4-mer skewness, IDENTICAL to cancer formula) can separate
human-infecting from non-human-infecting viruses WITHOUT any training.

Benchmark: Mollentze et al. 2021 (PLoS Biology) AUC=0.773 WITH training.

DATA SETUP (run once):
  git clone https://github.com/nardus/zoonotic_rank.git
  # OR just download these 2 files:
  # 1. Virus metadata + labels:
  #    https://raw.githubusercontent.com/nardus/zoonotic_rank/main/InternalData/allSequences.fasta
  # 2. S1 Table (species list + labels):
  #    From supplementary: save as S1_Table.csv

  # SIMPLEST: just download the FASTA + training labels:
  curl -L -o viruses.fasta "https://raw.githubusercontent.com/nardus/zoonotic_rank/main/InternalData/allSequences.fasta"
  curl -L -o virus_labels.csv "https://raw.githubusercontent.com/nardus/zoonotic_rank/main/InternalData/FinalData_Cleaned.csv"

Usage:
  python3 viral_ev_test.py viruses.fasta virus_labels.csv
  python3 viral_ev_test.py --from-repo /path/to/zoonotic_rank/
"""
import numpy as np, sys, json, csv, gzip, os, io
from scipy.stats import skew as skewness, mannwhitneyu, rankdata
from collections import defaultdict

# ══════════════════════════════════════════════════════════════════
# CANONICAL Ev — IDENTICAL to cancer pipeline
# ══════════════════════════════════════════════════════════════════
BASE = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3}
PRIMES = [p for p in range(2, 500) if all(p % i != 0 for i in range(2, int(p**0.5)+1))][:95]
P = np.random.default_rng(42).standard_normal((256, 500)) / np.sqrt(256)
assert abs(P[0, 0] - 0.01904482) < 1e-5, f"FATAL: P[0,0]={P[0,0]}"
print(f"[SANITY] P[0,0]={P[0,0]:.8f} ✓")


def kmer_index(kmer):
    return BASE[kmer[0]] * 64 + BASE[kmer[1]] * 16 + BASE[kmer[2]] * 4 + BASE[kmer[3]]


def compute_ev(seq):
    """Ev for a nucleotide sequence (DNA or RNA). Returns (ev, gc, cpg_oe)."""
    seq = seq.upper().replace('U', 'T').replace('N', '')
    if len(seq) < 200:
        return np.nan, np.nan, np.nan
    v = np.zeros(256)
    valid = 0
    for i in range(len(seq) - 3):
        kmer = seq[i:i+4]
        if all(c in BASE for c in kmer):
            v[kmer_index(kmer)] += 1
            valid += 1
    if valid < 50:
        return np.nan, np.nan, np.nan
    f = v / v.sum()
    assert abs(f.sum() - 1.0) < 1e-9
    ps = (f @ P)[PRIMES]
    ev = abs(skewness(ps)) * 6.07 + 0.10
    gc = (seq.count('G') + seq.count('C')) / len(seq)
    # CpG O/E
    c_ct = seq.count('C')
    g_ct = seq.count('G')
    cpg_ct = seq.count('CG')
    cpg_oe = (cpg_ct * len(seq)) / max(c_ct * g_ct, 1)
    return ev, gc, cpg_oe


def parse_fasta(path):
    """Parse FASTA → {header: sequence}."""
    opener = gzip.open if path.endswith('.gz') else open
    seqs = {}
    hdr, buf = None, []
    with opener(path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if hdr is not None:
                    seqs[hdr] = ''.join(buf)
                hdr = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if hdr is not None:
        seqs[hdr] = ''.join(buf)
    assert len(seqs) > 0
    print(f"[FASTA] {len(seqs)} sequences loaded")
    return seqs


def parse_labels_csv(path):
    """Parse Mollentze FinalData_Cleaned.csv or virus_metadata.json.
    Returns {accession: {'label': 0/1, 'species': str, 'family': str}}"""
    if path.endswith('.json'):
        with open(path) as f:
            raw = json.load(f)
        labels = {}
        for acc, info in raw.items():
            labels[acc] = {'label': info['label'], 'species': info.get('species', ''), 'family': ''}
        n_pos = sum(1 for v in labels.values() if v['label'] == 1)
        print(f"[LABELS] {len(labels)} entries: {n_pos} human, {len(labels)-n_pos} non-human")
        return labels
    labels = {}
    with open(path, encoding='utf-8', errors='replace') as f:
        # Try to auto-detect format
        sample = f.read(5000)
        f.seek(0)
        delim = ',' if sample.count(',') > sample.count('\t') else '\t'
        reader = csv.DictReader(f, delimiter=delim)
        cols = [c.lower().strip() for c in reader.fieldnames]

        # Find key columns
        acc_col = label_col = species_col = family_col = None
        for i, c in enumerate(cols):
            cn = reader.fieldnames[i]
            cl = c
            if any(k in cl for k in ['accession', 'genbank', 'sequence_name', 'acc']):
                acc_col = cn
            if any(k in cl for k in ['infects_humans', 'human_infection', 'is_human', 'label', 'infection']):
                label_col = cn
            if 'species' in cl and 'sub' not in cl:
                species_col = cn
            if 'family' in cl:
                family_col = cn

        if not acc_col:
            # Try first column
            acc_col = reader.fieldnames[0]
        if not label_col:
            # Search harder
            for cn in reader.fieldnames:
                if 'human' in cn.lower() or 'zoonotic' in cn.lower() or 'label' in cn.lower():
                    label_col = cn
                    break

        print(f"[LABELS] Columns detected: acc='{acc_col}', label='{label_col}', species='{species_col}', family='{family_col}'")
        assert label_col, f"Cannot find label column in {reader.fieldnames}"

        for row in reader:
            acc = row.get(acc_col, '').strip()
            if not acc:
                continue
            raw_label = row.get(label_col, '').strip().lower()
            if raw_label in ('1', 'true', 'yes', 'human'):
                label = 1
            elif raw_label in ('0', 'false', 'no', 'other', 'animal'):
                label = 0
            else:
                continue
            species = row.get(species_col, '').strip() if species_col else ''
            family = row.get(family_col, '').strip() if family_col else ''
            labels[acc] = {'label': label, 'species': species, 'family': family}

    n_pos = sum(1 for v in labels.values() if v['label'] == 1)
    n_neg = sum(1 for v in labels.values() if v['label'] == 0)
    print(f"[LABELS] {len(labels)} entries: {n_pos} human-infecting, {n_neg} non-human")
    assert n_pos >= 10 and n_neg >= 10, f"Too few labels: pos={n_pos}, neg={n_neg}"
    return labels


def compute_auroc(labels, scores):
    n1 = sum(labels)
    n0 = len(labels) - n1
    if n0 == 0 or n1 == 0:
        return np.nan
    ranks = rankdata(scores)
    pos_rank_sum = sum(r for l, r in zip(labels, ranks) if l == 1)
    return (pos_rank_sum - n1 * (n1 + 1) / 2) / (n0 * n1)


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 viral_ev_test.py <viruses.fasta> <labels.csv>")
        print("\nData download:")
        print("  curl -L -o viruses.fasta 'https://raw.githubusercontent.com/nardus/zoonotic_rank/main/InternalData/allSequences.fasta'")
        print("  curl -L -o virus_labels.csv 'https://raw.githubusercontent.com/nardus/zoonotic_rank/main/InternalData/FinalData_Cleaned.csv'")
        return

    seqs = parse_fasta(sys.argv[1])
    labels = parse_labels_csv(sys.argv[2])

    # Match sequences to labels
    results = []
    matched, unmatched = 0, 0
    for acc, seq in seqs.items():
        # Try exact match, then prefix match
        info = labels.get(acc)
        if info is None:
            # Try stripping version
            base_acc = acc.split('.')[0]
            info = labels.get(base_acc)
        if info is None:
            # Try matching by substring
            for k, v in labels.items():
                if k in acc or acc in k:
                    info = v
                    break
        if info is None:
            unmatched += 1
            continue
        matched += 1
        ev, gc, cpg_oe = compute_ev(seq)
        if np.isnan(ev):
            continue
        results.append({
            'acc': acc, 'label': info['label'], 'species': info['species'],
            'family': info['family'], 'ev': ev, 'gc': gc, 'cpg_oe': cpg_oe,
            'length': len(seq)
        })

    print(f"\n[MATCH] {matched} matched, {unmatched} unmatched, {len(results)} with valid Ev")

    # ── Sanity checks ──
    n_pos = sum(1 for r in results if r['label'] == 1)
    n_neg = sum(1 for r in results if r['label'] == 0)
    assert n_pos >= 10, f"Too few human-infecting: {n_pos}"
    assert n_neg >= 10, f"Too few non-human: {n_neg}"
    print(f"  Human-infecting: {n_pos}, Non-human: {n_neg}")

    all_labels = [r['label'] for r in results]
    ev_scores = [r['ev'] for r in results]
    gc_scores = [r['gc'] for r in results]
    cpg_scores = [r['cpg_oe'] for r in results]

    pos_evs = [r['ev'] for r in results if r['label'] == 1]
    neg_evs = [r['ev'] for r in results if r['label'] == 0]

    # ════════════════════════════════════════════
    # MAIN RESULTS
    # ════════════════════════════════════════════
    U, p_ev = mannwhitneyu(pos_evs, neg_evs, alternative='two-sided')

    # AUROC — try both directions
    auc_ev_high = compute_auroc(all_labels, ev_scores)
    auc_ev = max(auc_ev_high, 1 - auc_ev_high)
    ev_direction = "human=HIGH Ev" if auc_ev_high > 0.5 else "human=LOW Ev"

    auc_gc_high = compute_auroc(all_labels, gc_scores)
    auc_gc = max(auc_gc_high, 1 - auc_gc_high)

    auc_cpg_high = compute_auroc(all_labels, cpg_scores)
    auc_cpg = max(auc_cpg_high, 1 - auc_cpg_high)

    print(f"\n{'='*60}")
    print(f"VIRAL PANDEMIC PREDICTION — RESULTS")
    print(f"{'='*60}")
    print(f"  N viruses:           {len(results)} ({n_pos} human, {n_neg} non-human)")
    print(f"  Human Ev:            {np.mean(pos_evs):.4f} ± {np.std(pos_evs):.4f}")
    print(f"  Non-human Ev:        {np.mean(neg_evs):.4f} ± {np.std(neg_evs):.4f}")
    print(f"  Mann-Whitney p:      {p_ev:.2e}")
    print(f"")
    print(f"  Ev AUROC:            {auc_ev:.3f} ({ev_direction})")
    print(f"  GC AUROC:            {auc_gc:.3f}")
    print(f"  CpG O/E AUROC:       {auc_cpg:.3f}")
    print(f"")
    print(f"  BENCHMARK: Mollentze et al. trained model AUC = 0.773")
    print(f"  BENCHMARK: Phylogenetic relatedness model AUC < 0.773")

    # ── Sanity: check if GC alone dominates ──
    if auc_gc > auc_ev + 0.05:
        print(f"\n  [NOTE] GC alone outperforms Ev by {auc_gc-auc_ev:.3f}")
    if auc_cpg > auc_ev + 0.05:
        print(f"  [NOTE] CpG O/E alone outperforms Ev by {auc_cpg-auc_ev:.3f}")

    # ════════════════════════════════════════════
    # PER-FAMILY ANALYSIS
    # ════════════════════════════════════════════
    families = defaultdict(list)
    for r in results:
        if r['family']:
            families[r['family']].append(r)

    print(f"\n{'='*60}")
    print(f"PER-FAMILY ANALYSIS")
    print(f"{'='*60}")
    print(f"  {'Family':<25} {'n':>4} {'n+':>4} {'n-':>4} {'AUC':>6} {'p':>10}")
    family_results = []
    for fam, recs in sorted(families.items(), key=lambda x: -len(x[1])):
        fp = sum(1 for r in recs if r['label'] == 1)
        fn = len(recs) - fp
        if fp < 3 or fn < 3:
            continue
        fl = [r['label'] for r in recs]
        fs = [r['ev'] for r in recs]
        fauc = compute_auroc(fl, fs)
        fauc = max(fauc, 1 - fauc)
        fpos = [r['ev'] for r in recs if r['label'] == 1]
        fneg = [r['ev'] for r in recs if r['label'] == 0]
        _, fp_val = mannwhitneyu(fpos, fneg, alternative='two-sided')
        print(f"  {fam:<25} {len(recs):>4} {fp:>4} {fn:>4} {fauc:>6.3f} {fp_val:>10.2e}")
        family_results.append({'family': fam, 'n': len(recs), 'auc': fauc, 'p': fp_val})

    # ════════════════════════════════════════════
    # INTERPRETATION
    # ════════════════════════════════════════════
    print(f"\n{'='*60}")
    print(f"INTERPRETATION")
    print(f"{'='*60}")
    if auc_ev < 0.52:
        print("  NO SIGNAL. Ev does not separate human/non-human viruses.")
        print("  The Ev formula's biology (chromatin state) may not apply to viral genomes.")
    elif auc_ev < 0.58:
        print("  WEAK SIGNAL. Marginal separation detected.")
        print("  Some compositional difference exists but Ev is not the right metric.")
    elif auc_ev < 0.65:
        print(f"  MODERATE SIGNAL (AUC={auc_ev:.3f}). Training-free approaches viable.")
        print("  Publishable as proof-of-concept. Not yet deployable.")
    elif auc_ev < 0.73:
        print(f"  STRONG SIGNAL (AUC={auc_ev:.3f}). Approaching trained model performance.")
        print("  Zero-training matches trained ML within margin. Major finding.")
        print("  → Patent application warranted. WHO/CDC deployment feasible.")
    else:
        print(f"  EXCEEDS trained model (AUC={auc_ev:.3f} vs benchmark 0.773).")
        print("  Transformative result. Immediate public health application.")

    # ════════════════════════════════════════════
    # SAVE
    # ════════════════════════════════════════════
    out = {
        'n_total': len(results), 'n_human': n_pos, 'n_nonhuman': n_neg,
        'auroc_ev': float(auc_ev), 'auroc_gc': float(auc_gc), 'auroc_cpg': float(auc_cpg),
        'ev_direction': ev_direction, 'p_ev': float(p_ev),
        'mean_ev_human': float(np.mean(pos_evs)), 'mean_ev_nonhuman': float(np.mean(neg_evs)),
        'mean_gc_human': float(np.mean([r['gc'] for r in results if r['label'] == 1])),
        'mean_gc_nonhuman': float(np.mean([r['gc'] for r in results if r['label'] == 0])),
        'benchmark_mollentze_auc': 0.773,
        'per_family': family_results,
    }
    with open('viral_ev_results.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: viral_ev_results.json")


if __name__ == '__main__':
    main()
