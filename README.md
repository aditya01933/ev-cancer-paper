# Ev Cancer Paper — Verified Project Files

## Paper
"A single formula predicts where cancer mutations and epigenetic aging converge from DNA sequence alone"
Aditya Tiwari, 2026
bioRxiv MS ID: BIORXIV/2026/710061

## Directory Structure
```
data/               Verified result JSONs (all cross-checked against manuscript)
raw_zones/          23 chromosome zone assignments + mutation-annotated versions
figures/            Publication figures (Fig1-6 from pipeline, Fig7 from Horvath analysis)
scripts/            Verified analysis scripts
manuscript/         Submitted documents
supplementary/      Additional File 1
validation/         Negative results and validation experiments
```

## Critical Data Files

| File | Contains | Verified |
|------|----------|----------|
| simpson_per_chr_result.json | 23 chr mutations + windows | OR=1.682 ✅ |
| pan_cancer_results.json | 15 cancer type ORs | mean=1.63 ✅ |
| two_mechanisms_corrected.json | ChIP M1=1.775x M2=2.287x | p=1.6e-39, 3.3e-99 ✅ |
| horvath_ev_zones_results.json | 346 clock CpGs mapped | OR=1.94 ✅ |
| cancer_methylation_results.json | 32 tumor-normal pairs | Δβ=-0.054 ✅ |
| germline_baseline_chr1_5.json | 1000G 2504 individuals | Z3/Z1=1.129 ✅ |
| zone3_survival_result.json | 997 BRCA patients | p=0.404 (negative) ✅ |
| drug_targets_zone3.json | 101 COSMIC genes | TP53=Z3, KRAS=Z1 ✅ |
| maths_rerun_results.json | Seed stability | 12/20 ✅ |
| cosmic_signature_results.json | Genome-wide 5-class | C>T OR=1.115 ✅ |

## WARNING: Do NOT Use These Files

| File | Location | Problem |
|------|----------|---------|
| two_mechanism_results.json | imp-research/ | Tests WRONG zone+mark combinations (Z1+H3K9me3 instead of Z3+H3K9me3) |
| exp_a3_per_chr_fixed.png | imp-research/ | Pre-reboot buggy session, inverted zones |
| exp1_per_cancer_or.png | imp-research/ | Broken rendering |
| exp5_germline_vs_somatic.png | imp-research/ | Shows Somatic OR=0.025 (pipeline error) |

## Ev Formula
```python
P = np.random.default_rng(42).standard_normal((256,500)) / sqrt(256)
# P[0,0] = 0.01904482
# 4-mer freq → P projection → skewness of 95 prime indices → scale
# Ev_resid = (Ev - (-5.938 * GC + 4.471)) / 0.456
# Zone 1: Ev_resid >= +0.382 | Zone 3: Ev_resid <= -1.471
```
