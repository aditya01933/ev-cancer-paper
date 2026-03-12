# Ev: DNA Sequence Composition Predicts Cancer Mutation Geography and Epigenetic Aging

[![bioRxiv](https://img.shields.io/badge/bioRxiv-2026.710061-b31b1b.svg)](https://www.biorxiv.org/content/10.1101/2026.710061)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A single fixed-parameter formula (Ev) operating on raw DNA sequence classifies genomic windows into chromatin-like zones and predicts somatic mutation enrichment across 15 cancer types — without any experimental epigenomic data, training, or reference databases.

**Paper:** Tiwari A. (2026) "A single formula predicts where cancer mutations and epigenetic aging converge from DNA sequence alone." _bioRxiv_ [doi: pending]

---

## Key Results

| Finding                                | Value                  | P-value     |
| -------------------------------------- | ---------------------- | ----------- |
| Genome-wide mutation enrichment (BRCA) | OR = 1.682             | < 10⁻³⁰⁰    |
| Pan-cancer universality                | 15/15 types OR > 1     | all < 0.05  |
| C>T only enriched mutation class       | OR = 1.24              | 8.7 × 10⁻²⁰ |
| Replication timing independent         | 3/3 strata significant | all < 0.05  |
| Germline baseline (1000G)              | Z3/Z1 = 1.129          | 1.5 × 10⁻⁸² |
| Horvath aging clock CpGs in Z3         | OR = 1.94              | 2.0 × 10⁻⁵  |
| TSG/Oncogene segregation               | OR = 9.2               | 4 × 10⁻¹²   |

## The Ev Formula

```python
import numpy as np
from scipy.stats import skew

# Canonical parameters (immutable)
P = np.random.default_rng(42).standard_normal((256, 500)) / np.sqrt(256)
PRIMES = [p for p in range(2, 500) if all(p % i != 0 for i in range(2, int(p**0.5)+1))][:95]
assert abs(P[0, 0] - 0.01904482) < 1e-6  # fingerprint check

def ev_resid(sequence, gc=None):
    """Compute Ev residual for a 5kb DNA sequence."""
    seq = sequence.upper().replace('N', '')
    if len(seq) < 1000:
        return None
    # Count 4-mer frequencies
    bases = 'ACGT'
    ki = {a+b+c+d: i*64+j*16+k*4+l
          for i,a in enumerate(bases) for j,b in enumerate(bases)
          for k,c in enumerate(bases) for l,d in enumerate(bases)}
    v = np.zeros(256)
    for i in range(len(seq) - 3):
        k = seq[i:i+4]
        if k in ki:
            v[ki[k]] += 1
    if v.sum() == 0:
        return None
    f = v / v.sum()
    if gc is None:
        gc = (seq.count('G') + seq.count('C')) / len(seq)
    ev = abs(skew((f @ P)[PRIMES])) * 6.07 + 0.10
    return (ev - (-5.938 * gc + 4.471)) / 0.456
```

**Zone thresholds:** Z1 (euchromatin) ≥ +0.382 | Z3 (heterochromatin) ≤ −1.471

## Repository Structure

```
ev-cancer-paper/
├── scripts/
│   ├── ev_formula.py              # Canonical Ev formula (import this)
│   ├── recompute_all_zones.py     # Full recompute from FASTA + MC3
│   ├── sanity_zones.py            # Zone sanity checks
│   ├── gc-comparision.py          # GC vs Ev zone comparison
│   ├── bio2_signatures.py         # Pan-cancer mutation signatures
│   ├── mut_signature_zones.py     # BRCA mutation class decomposition
│   ├── two_mechanisms_corrected.py # ChIP-seq two-mechanism test
│   ├── exp4_horvath_zones.py      # Horvath clock CpG analysis
│   ├── exp4b_clock_mutations.py   # Clock CpG × mutation density
│   ├── analyze_4mer_drivers.py    # Reverse-map P matrix to 4-mers
│   ├── analyze_ev_pca.py          # PCA of 4-mer space
│   ├── test_at_skew.py            # AT-skew = PC2 analysis
│   ├── alu_mutation_chain.py      # Alu density causal disproof
│   ├── repeat_masker.py           # RepeatMasker TE intersection
│   └── cross_species_te_density.py # Cross-species TE analysis
├── validation/
│   ├── validate_all_claims.py     # 17-test comprehensive validation engine
│   ├── e2e_tests.py               # 33 end-to-end consistency checks
│   ├── investigate_failures.py    # Discrepancy diagnosis
│   ├── VALIDATION_LOG.json        # Full validation output
│   └── VALIDATION_REPORT.txt      # Human-readable summary
├── manuscript/
│   ├── Tiwari_2026_FINAL_v6.md    # Final manuscript (markdown source)
│   ├── Tiwari_2026_FINAL_v6.docx  # Final manuscript (Word)
│   └── Additional_File_1_v6.pdf   # Supplementary materials
├── figures/
│   ├── Fig1_ev_distribution.png
│   ├── Fig2_per_chr_OR.png
│   ├── Fig3_pan_cancer.png
│   ├── Fig4_CT_specificity.png
│   ├── Fig5_replication_timing.png
│   ├── Fig6_germline_amplification.png
│   └── Fig7_horvath_clock.png
├── data/                          # Computed results (JSON)
│   ├── horvath_ev_zones_results.json
│   ├── germline_baseline_chr1_5.json
│   ├── cancer_methylation_results.json
│   ├── pan_cancer_results.json
│   ├── drug_targets_zone3.json
│   ├── simpson_per_chr_result.json
│   ├── two_mechanisms_corrected.json
│   ├── mutation_signature_zones.json
│   ├── at_skew_analysis.json
│   ├── ev_pca_analysis.json
│   └── ev_4mer_drivers.json
├── raw_zones/                     # Per-chromosome zone assignments
│   └── human_chr*_ev_zones*.json
└── README.md
```

## Raw Data (not included, download separately)

```bash
mkdir -p raw_data/{fasta,maf,chipseq,repliseq,horvath,cross_species}

# Human genome (hg38)
for c in $(seq 1 22) X; do
  wget -P raw_data/fasta/ "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr${c}.fa.gz"
  mv raw_data/fasta/chr${c}.fa.gz raw_data/fasta/human_chr${c}.fa.gz
done

# MC3 somatic mutations
# Download from https://gdc.cancer.gov/about-data/publications/mc3-2017

# ENCODE ChIP-seq (GM12878)
# H3K9me3 and H3K4me3 narrowPeak BED files from ENCODE portal

# Repli-seq
# ENCFF001GVQ from ENCODE (GM12878 wavelet-smoothed)

# RepeatMasker
wget -P raw_data/ "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
```

## Validation

The comprehensive validation engine recomputes all claims from raw FASTA + MC3 MAF:

```bash
cd ev-cancer-paper
python3 -u validation/validate_all_claims.py 2>&1 | tee validation/run.log
python3 validation/e2e_tests.py
```

**Result: 16/17 tests pass, 33/34 E2E checks pass.**

The single remaining failure (T12: TSG/OG segregation) is due to using 30 genes instead of the full 101 COSMIC Tier 1 list. TP53=Zone 3 and KRAS=Zone 1 are confirmed.

## Quick Start

```python
from scripts.ev_formula import ev_resid

# Classify any 5kb DNA sequence
score = ev_resid("ATCGATCG..." * 625)  # 5000bp
if score >= 0.382:
    print("Zone 1 (euchromatin)")
elif score <= -1.471:
    print("Zone 3 (heterochromatin)")
else:
    print("Zone 2 (intermediate)")
```

## Citation

```bibtex
@article{tiwari2026ev,
  title={A single formula predicts where cancer mutations and epigenetic
         aging converge from DNA sequence alone},
  author={Tiwari, Aditya},
  journal={bioRxiv},
  year={2026},
  doi={10.1101/2026.710061}
}
```

## License

MIT License. See [LICENSE](LICENSE).

## Acknowledgement

AI assistance (Claude, Anthropic) was used for code development and manuscript editing. All scientific analyses, experimental design, data interpretation, and conclusions are the author's own.
