# A single formula predicts where cancer mutations and epigenetic aging converge from DNA sequence alone

Aditya Tiwari^1\*^

^1^ Independent Researcher, Goa, India

\*Correspondence: adityatiwari01933@gmail.com

## Abstract

Regional variation in somatic mutation rates across cancer genomes is strongly associated with chromatin organization, but existing predictive models require experimental epigenomic data such as ChIP-seq or Hi-C. Here I show that a single formula operating on raw DNA sequence---a 4-mer frequency skewness statistic termed Ev---classifies genomic windows into chromatin-like zones without any experimental input, training data, or reference databases. Applying Ev to 582,028 non-overlapping 5-kb windows across the human genome, I find that the lowest-scoring zone (Zone 3, corresponding to heterochromatin) carries 1.68-fold more somatic mutations than the highest-scoring zone (Zone 1, euchromatin) in 992 TCGA breast cancer genomes (odds ratio = 1.682, P < 10^-300^). This enrichment is universal across all 15 TCGA cancer types tested (mean OR = 1.63, range 1.43--1.91, zero exceptions) and persists after stratification by replication timing (P < 0.05 in all strata). The effect is driven specifically by C>T transitions consistent with SBS1 (5-methylcytosine deamination; OR = 1.24, P = 8.7 x 10^-20^), and is independently confirmed by H3K4me3 ChIP-seq enrichment analysis (3.45-fold, P = 7.2 x 10^-223^). Zone 3 vulnerability is constitutional: it is present in healthy germline genomes from 2,504 individuals (Z3/Z1 = 1.13, P = 1.5 x 10^-82^) and amplified 1.38-fold in cancer. Sequence composition alone also discriminates tumor suppressors from oncogenes among 101 COSMIC Tier 1 cancer driver genes (OR = 9.2, P = 4 x 10^-12^). Moreover, Zone 3 is enriched 1.94-fold (P = 2.0 x 10^-5^) for the CpG sites that compose Horvath's multi-tissue epigenetic aging clock, specifically those that gain methylation with age, linking cancer mutation geography and epigenetic aging through shared sequence-encoded chromatin vulnerability. These results demonstrate that the regional distribution of somatic mutations across cancer genomes can be predicted from DNA sequence composition alone, without recourse to experimental epigenomics.

**Keywords:** somatic mutation rate, chromatin, sequence composition, 4-mer, heterochromatin, cancer genomics, TCGA, epigenetic aging

## Background

Somatic mutations in cancer genomes are not distributed uniformly. Instead, mutation rates vary by more than five-fold across megabase-scale regions within individual tumor genomes [1--3]. Understanding the determinants of this variation has important implications for distinguishing driver from passenger mutations, for interpreting non-coding variation, and for modeling cancer evolution.

Chromatin organization has emerged as a dominant factor explaining this regional variation. In a landmark study, Schuster-Bockler and Lehner showed that levels of the heterochromatin-associated histone modification H3K9me3 alone can account for more than 40% of mutation rate variation at the megabase scale [1]. Subsequent work has extended these observations to include replication timing [4], transcription factor binding [5], nucleosome positioning [6], and three-dimensional genome topology [7]. Integrative models combining multiple epigenomic features can explain over 55% of regional mutation rate variance [1, 8].

A common thread in this literature is the requirement for experimental epigenomic data. Predicting where mutations will accumulate in a given genome requires, at minimum, ChIP-seq profiles, replication timing assays, Hi-C contact maps, or curated databases of regulatory features. A recent study by Liu et al. used DNA motifs as surrogates for epigenetic signals to predict somatic mutation rates across 13 cancers, but still required a trained neural network model and transcription factor motif databases [9]. Whether regional mutation rate variation can be predicted from DNA sequence composition alone---without any experimental input---has remained an open question. Liu et al. noted explicitly that they were "not aware of any study that can predict mutation rates at even kilobase resolution using only sequences" [9].

Here I address this question. I describe a fixed-parameter formula ("Ev") that operates on 4-mer frequencies in raw DNA sequence to classify 5-kb genomic windows into zones corresponding to euchromatin and heterochromatin. I then demonstrate that the heterochromatin-associated zone (Zone 3) is universally enriched for somatic mutations across 15 cancer types, and I characterize the mutational mechanism and its constitutional basis.

## Results

### A sequence-only formula classifies genomic windows into chromatin zones

The Ev statistic is computed from the frequencies of all 256 tetranucleotides (4-mers) in a 5-kb genomic window. These frequencies are projected through a fixed random matrix (256 x 500, generated from a standard normal distribution with seed 42), and the absolute skewness of the 95 prime-indexed projections, linearly scaled, yields the Ev score (see Methods). To remove the dominant effect of GC content, I regress Ev against GC fraction across all windows and use the residual (Ev_resid) for zone classification.

Three zones are defined by fixed thresholds on Ev_resid: Zone 1 (Ev_resid >= +0.382; euchromatin-like, GC-rich, gene-dense), Zone 2 (intermediate), and Zone 3 (Ev_resid <= -1.471; heterochromatin-like, AT-rich, transposable element-dense). Across the human genome (chromosomes 1--22 and X), this partitioning yields 152,519 Zone 1 windows (26.2%), 287,733 Zone 2 windows (49.4%), and 141,776 Zone 3 windows (24.4%) out of 582,028 total.

It is important to note what Ev is and what it is not. By principal component analysis, the first component of 4-mer frequency space in human genomic windows explains 76.5% of variance and aligns with the AT/GC composition axis. The correlation between Ev and GC content is r = 0.625, and even after residualization, Ev's area under the receiver operating characteristic curve (AUC) for binary chromatin-state prediction is 0.647, compared to 0.850 for GC content alone and 0.874 for an optimized linear discriminant. Thus Ev is a noisy, partially linear proxy for GC-related sequence features. Its utility lies not in mathematical novelty but in its capacity to extract a biologically meaningful chromatin signal from raw sequence, as I demonstrate below.

### Zone 3 predicts somatic mutation enrichment genome-wide in breast cancer

I obtained open-access somatic mutation calls (masked MAF files) for 992 TCGA-BRCA patients from the Genomic Data Commons (GDC). After intersecting 89,565 mutations with the zone map, mutation rates per window were: Zone 1, 0.121; Zone 2, 0.155; Zone 3, 0.188 (Table 1). The ratio of Zone 3 to Zone 1 mutation rates was 1.554, with a Fisher exact test odds ratio of 1.682 (P < 10^-300^, one-sided). Zone 2 fell between Zone 1 and Zone 3, confirming a monotonic dose-response relationship.

This result was reproducible at the single-chromosome level: chromosome 1 alone gave OR = 1.680 (P = 3.6 x 10^-59^), nearly identical to the genome-wide estimate. To assess consistency across all chromosomes, I computed per-chromosome odds ratios (Fig. 2). Of 23 chromosomes tested, 22 showed OR > 1, and 21 reached statistical significance at P < 0.05 (median OR = 1.44, range 0.90--2.24). The sole exception was chromosome 19, where Zone 3 showed a marginally inverted OR of 0.90. Chromosome 19 is the most gene-dense human chromosome, with GC content of ~48% compared to a genome average of 41%; consequently, most of chromosome 19 classifies as Zone 3 by sequence composition even though it is functionally euchromatic. The highest per-chromosome OR was observed on chromosome X (OR = 2.24, P = 2.2 x 10^-71^), consistent with X-inactivation creating large blocks of constitutive heterochromatin.

### Pan-cancer universality across 15 TCGA cancer types

To test whether Zone 3 enrichment is specific to breast cancer or reflects a universal property of somatic mutagenesis, I repeated the analysis across 15 TCGA cancer types (50 randomly sampled patients per type; see Methods). Zone 3 mutation enrichment was significant in every cancer type tested (Table 2; Fig. 3). Odds ratios ranged from 1.43 (GBM, glioblastoma) to 1.91 (PRAD, prostate adenocarcinoma), with a mean of 1.63 and zero exceptions to the Z3 > Z2 > Z1 dose-response gradient. This universality held across cancers dominated by different mutational processes: tobacco-associated lung adenocarcinoma (LUAD, OR = 1.52), UV-associated melanoma (SKCM, OR = 1.77), mismatch repair-deficient uterine cancer (UCEC, OR = 1.80, P < 10^-300^), and low-mutation-rate thyroid cancer (THCA, OR = 1.57, P = 4.5 x 10^-5^).

### C>T transitions are the only mutation class enriched in Zone 3

To identify which mutational process drives Zone 3 enrichment, I decomposed BRCA somatic mutations by substitution class and computed zone-specific rates for each (Table 3). C>T transitions were the only class with a Z3/Z1 ratio exceeding 1.0 (ratio = 1.12, OR = 1.24, P = 8.7 x 10^-20^). All other classes were either neutral (C>A, C>G) or depleted in Zone 3 relative to Zone 1: T>A (OR = 0.78, P = 7.9 x 10^-6^), T>C (OR = 0.74, P = 1.8 x 10^-11^), and T>G (OR = 0.70, P = 7.1 x 10^-9^).

The specific enrichment of C>T in heterochromatin-associated zones is consistent with COSMIC mutational signature SBS1, which reflects spontaneous deamination of 5-methylcytosine (5mC) at CpG dinucleotides [10]. CpG-dense regions in heterochromatin are constitutively methylated; 5mC undergoes spontaneous hydrolytic deamination to thymine, producing C>T transitions that are poorly repaired in condensed chromatin [11, 12]. The depletion of T-origin mutations (T>A, T>C, T>G) in Zone 3 is consistent with the AT-rich composition of heterochromatin, where fewer T bases are available per window once the dominant A-tract sequences are accounted for, and with more efficient repair of mismatches in transcriptionally active euchromatin [13].

This pattern was not restricted to breast cancer. In a pan-cancer signature analysis across LUAD (lung), SKCM (melanoma), and PRAD (prostate), C>T was the only class enriched in Zone 3 in all three additional cancer types, confirming C>T as the universal mechanism underlying Zone 3 mutation enrichment [see Additional file 1].

### ChIP-seq validation connects sequence-defined zones to chromatin biology

To connect the sequence-defined zones to established chromatin biology, I intersected Zone 1 and Zone 3 windows with ChIP-seq peaks for H3K9me3 (a heterochromatin mark) and H3K4me3 (an active promoter mark) from the GM12878 lymphoblastoid cell line (ENCODE).

Within Zone 3 (heterochromatin by sequence), windows overlapping H3K4me3 peaks carried 3.45-fold more mutations than non-overlapping windows (P = 7.2 x 10^-223^, chromosome 1). This striking enrichment confirms that even within the heterochromatin-associated compartment, the presence of active transcription marks---indicative of isolated actively transcribed loci embedded in otherwise condensed chromatin---further elevates mutation rates, consistent with transcription-coupled mutagenesis in which transient single-stranded DNA exposure increases vulnerability to DNA damage [14, 15].

The net mutation enrichment in Zone 3 (OR = 1.68 genome-wide) arises because repair failure in heterochromatin---which operates specifically through C>T accumulation---produces a larger cumulative signal than localized transcription-coupled damage in euchromatin, which is distributed across multiple mutation classes.

A caveat applies: the H3K4me3 ChIP-seq data covered chromosome 1 only, while the H3K9me3 data were genome-wide. Additionally, the GM12878 cell line is lymphoblastoid, not mammary. Cell-type matched ChIP-seq (e.g., MCF-7 H3K9me3) would further strengthen these observations.

### Zone 3 enrichment is independent of replication timing

Late replication is a known predictor of elevated mutation rates in cancer genomes [4, 16], and late-replicating regions overlap substantially with heterochromatin. To test whether Zone 3 enrichment is merely a proxy for replication timing, I obtained GM12878 Repli-seq data from ENCODE (ENCFF001GVQ) and stratified chromosome 1 windows into replication timing tertiles (early, mid, late) before computing zone-specific mutation rates within each tertile (Table 4).

Zone 3 carried significantly more mutations than Zone 1 in all three replication timing strata: early-replicating (OR = 1.37, P = 5.5 x 10^-5^), mid-replicating (OR = 1.28, P = 6.8 x 10^-4^), and late-replicating (OR = 1.16, P = 3.8 x 10^-2^). The cross-tabulation of zones against replication timing tertiles showed approximately uniform distribution: Zone 3 windows were not concentrated in late-replicating regions but spread evenly across all three timing classes (~3,200--3,300 windows each). The Ev-derived zone effect is therefore independent of replication timing.

Notably, the Ev zone effect was strongest in early-replicating regions (OR = 1.37) and weakest in late-replicating regions (OR = 1.16). This gradient is expected: in late-replicating windows, replication-associated errors already elevate mutation rates broadly, partially masking the heterochromatin-specific signal. In early-replicating windows, replication timing contributes less to the overall mutation rate, allowing the sequence-composition signal to dominate.

### Zone 3 vulnerability is constitutional and amplified in cancer

If Zone 3 enrichment reflects a physical property of heterochromatin rather than a cancer-specific process, it should also be detectable in healthy germline genomes. Using high-coverage variant calls from the 1000 Genomes Project (2,504 individuals, chromosomes 1--5, 26.2 million variants), I found that germline variant density was significantly higher in Zone 3 than Zone 1 (Z3/Z1 = 1.129, P = 1.5 x 10^-82^), with a perfect dose-response: Zone 1 (120.7 variants/window) < Zone 2 (125.1) < Zone 3 (136.2).

Comparing the germline Z3/Z1 ratio (1.13) to the somatic cancer ratio (1.55), cancer amplifies the constitutional Zone 3 vulnerability by a factor of 1.38. This amplification is consistent with a model in which methylation loss in tumor heterochromatin [17] exposes additional 5mC sites to deamination, accelerating the baseline C>T accumulation rate.

Supporting this model, analysis of matched tumor-normal methylation arrays (32 TCGA-BRCA pairs, Illumina 450K) showed that Zone 3 windows undergo significant hypomethylation in tumors (mean Delta-beta = -0.054, P = 3.5 x 10^-9^, observed in 91% of tumors), while Zone 1 windows show modest hypermethylation (Delta-beta = +0.032, P = 5.5 x 10^-7^). This is consistent with the classic cancer epigenome pattern of global heterochromatin hypomethylation coupled with focal euchromatin hypermethylation [18].

### Cross-species conservation of Zone 3

If the Ev formula captures a fundamental property of DNA sequence organization, it should detect similar zones in distantly related species. Scanning chromosome-scale sequences from ten eukaryotic species spanning approximately 1.5 billion years of evolution---including vertebrates (human, chimpanzee), dicots (Arabidopsis, soybean), grasses (maize, sorghum, Brachypodium), a moss (Physcomitrium), a lycophyte (Selaginella), and a chlorophyte alga (Chlamydomonas)---I found that Zone 3 was detected in all ten species with consistent Ev_resid characteristics. Across the nine land plant and vertebrate species (excluding the unicellular alga Chlamydomonas, which has minimal constitutive heterochromatin), the mean Zone 3 Ev_resid value was -2.44 with a coefficient of variation of 4.2%.

In Arabidopsis, where dense functional annotations are available, Zone 3 showed strong overlap with transposable element-dense pericentromeric regions (OR = 3.44 for windows with >59% TE content, P < 10^-50^) and reduced CpG observed-to-expected ratios (Zone 1 CpG O/E = 0.77 vs Zone 3 = 0.61), consistent with decades of CpG depletion through 5mC deamination in methylated heterochromatin.

An intriguing feature of the cross-species analysis is that the AT/GC asymmetry defining Zone 3 inverts between genomic kingdoms. In vertebrates and dicots, Zone 3 is AT-rich (driven by alpha-satellite repeats and Gypsy LTR transposable elements, respectively). In grass monocots and the high-GC alga Chlamydomonas, Zone 3 is relatively GC-rich against an AT-rich background (driven by Helitron transposable elements). Despite this polarity flip, the Ev formula detects Zone 3 in all cases because it captures the 4-mer skewness away from the genomic mean, regardless of which direction the skewness takes.

### Zone 3 enrichment extends to epigenetic aging clock CpGs

The results above establish that Zone 3 is a constitutionally mutation-prone compartment. To test whether the same sequence-defined zones are relevant to a different biological process---epigenetic aging---I examined the distribution of the 353 CpG sites that compose Horvath's multi-tissue DNA methylation age predictor [20]. These clock CpGs were selected by elastic net regression across 8,000 samples from 51 tissue types to predict chronological age from methylation levels, with no reference to chromatin state or sequence composition [20]. Of the 353 clock CpGs, 193 gain methylation with age (positive markers, enriched at Polycomb group targets) and 160 lose methylation with age (negative markers, enriched at CpG shores) [20].

Mapping the 353 clock CpGs to the Ev zone framework (346 successfully mapped; 7 fell in unmappable regions), I found a highly non-random distribution across zones (chi-squared = 20.77, P = 3.1 x 10^-5^). Clock CpGs were significantly depleted from Zone 1 (15.9% vs 26.2% genome-wide; Fisher OR = 0.53, P = 5.5 x 10^-6^) and significantly enriched in Zone 3 (30.6% vs 24.4%; Fisher OR = 1.37, P = 8.4 x 10^-3^).

Strikingly, the zone enrichment was driven almost entirely by the positive (age-hypermethylated) clock CpGs. These 190 successfully mapped CpGs showed 1.58-fold enrichment in Zone 3 (38.4% vs 24.4% genome-wide; Fisher OR = 1.94, P = 2.0 x 10^-5^), with corresponding depletion from Zone 1 (16.8% vs 26.2%; OR = 0.57, P = 2.9 x 10^-3^). By contrast, the 156 negative (age-hypomethylated) CpGs showed no significant Zone 3 enrichment (21.2% vs 24.4%; OR = 0.83, P = 0.40) and instead concentrated in Zone 2 (64.1% vs 49.4%).

This separation is mechanistically coherent. The positive clock CpGs---those that gain methylation with age---are known to reside near Polycomb group protein targets [20]. Polycomb-mediated repression maintains these CpGs in an unmethylated state during youth; as epigenetic maintenance erodes with age, the underlying sequence context (Zone 3, heterochromatin-associated, methylation-prone) drives progressive methylation gain. The negative clock CpGs, which lose methylation, sit in intermediate compositional zones where active maintenance is required to preserve methylation, and its age-related decline leads to methylation loss.

These findings connect the cancer mutation landscape to epigenetic aging through a shared sequence-encoded vulnerability: Zone 3 windows are both the sites where somatic mutations preferentially accumulate (OR = 1.682) and the sites where the epigenetic clock preferentially ticks (OR = 1.94 for age-gaining CpGs). Both processes reflect the consequences of impaired maintenance in heterochromatin---mutations from unrepaired 5mC deamination, and age-associated methylation gain from Polycomb erosion---and both are predictable from DNA sequence composition alone.

### Drug target genes segregate by zone

Classification of 101 COSMIC Cancer Gene Census Tier 1 genes by their majority zone assignment revealed a biologically coherent pattern. Tumor suppressors that are inactivated by accumulated mutations in cancer clustered in Zone 3: TP53 (100% Zone 3 windows), BRCA1 (72.7% Zone 3), ERBB2 (66.7%), ARID1A, STAT3, MCL1, and MYCN. By contrast, oncogenes activated by specific point mutations clustered in Zone 1: KRAS, NRAS, and HRAS (all three RAS family oncogenes).

This separation is consistent with the functional logic of cancer driver genes. Tumor suppressors require only loss of function, which can occur through the accumulation of any disruptive mutation---a process favored in the mutationally permissive Zone 3 environment. Oncogenes require specific activating point mutations at defined hotspot residues, a process independent of regional mutation rate.

Quantifying this separation across all 101 COSMIC Tier 1 genes, 73% of tumor suppressors fell in Zone 3 compared to 22% of oncogenes (Fisher OR = 9.2, P = 4 x 10^-12^). This effect size is 5.5-fold larger than the genome-wide mutation enrichment (OR = 1.68), indicating that sequence composition not only predicts regional mutation rates but also discriminates functional classes of cancer driver genes without labels or training data.

### Negative results

Not all tested associations were positive. Zone 3 burden did not predict overall survival in 997 TCGA-BRCA patients (log-rank P = 0.404). The Ev formula is therefore not a prognostic biomarker. Patients with high tumor mutational burden (TMB) showed lower Zone 3 enrichment, indicating that hypermutation dilutes the Zone 3 signal by introducing mutations in non-heterochromatic regions. Additionally, the Ev statistic does not approximate Kolmogorov complexity (correlation with LZMA compression ratio: r = 0.070).

## Discussion

These results establish that the regional distribution of somatic mutations across cancer genomes can be predicted from DNA sequence composition alone. A fixed-parameter formula, with no training and no experimental input, identifies genomic windows enriched for somatic mutations at OR = 1.68 (P < 10^-300^) in breast cancer and with complete universality across 15 TCGA cancer types.

The finding itself---that heterochromatin accumulates more mutations than euchromatin---is not new. Schuster-Bockler and Lehner demonstrated this relationship in 2012 using H3K9me3 ChIP-seq [1], and it has been confirmed and extended by many subsequent studies [3, 8, 9, 16]. What is new is the input: raw FASTA sequence, nothing else. The entire analytical pipeline requires only a reference genome assembly and somatic mutation coordinates. No ChIP-seq, no Hi-C, no methylation arrays, no trained classifier, no curated feature database.

A natural objection is that Ev is largely a proxy for GC content (r = 0.625), and GC content itself correlates with chromatin state [19]. This is correct. I have shown that GC content alone achieves higher AUC (0.850) than Ev (0.647) for binary chromatin classification. The contribution of this work is not that Ev outperforms GC content---it does not---but that even a crude sequence composition measure, applied without optimization, produces a robust and universal cancer mutation signal. The signal arises because the underlying biology is strong: constitutive heterochromatin methylation drives CpG deamination at rates high enough to produce genome-wide enrichment detectable across hundreds of patients. The Ev formula, for all its noise, captures enough of this biology to produce OR = 1.68 from sequence alone.

Reverse-mapping the random projection matrix reveals which sequence features drive the Ev signal. Analytically, the 256 x 500 projection matrix assigns approximately equal weight to all 4-mers (CpG vs non-CpG L2 norm difference: <2%), confirming that the biological signal arises from input frequency distributions, not projection bias. Empirically, Ev_resid correlates most strongly with poly-A tract frequency (partial r = +0.41 for AAAT after GC control) and anti-correlates with T-rich runs (partial r = -0.40 for CTTT). CpG-containing 4-mers are less correlated with Ev_resid than non-CpG 4-mers (mean |r| = 0.167 vs 0.203), indicating that the formula does not preferentially detect CpG content. A linear model using the top 20 driving 4-mers achieves R-squared = 0.40 against Ev_resid, confirming that the random projection captures nonlinear interactions among 4-mer frequencies that cannot be replaced by a simple deterministic formula.

Principal component analysis of the 256-dimensional 4-mer frequency space clarifies this further. The first principal component (35.3% of variance) aligns almost perfectly with GC content (r = -0.994) and has zero correlation with Ev_resid (r = +0.003). Ev_resid instead correlates with the second principal component (16.9% of variance, r = +0.380), which represents strand compositional asymmetry: AT-skew, defined as (A-T)/(A+T). Zone 1 windows are A-enriched on the reference strand (AT-skew = +0.041), while Zone 3 windows are T-enriched (AT-skew = -0.048), consistent with transcription-coupled repair creating strand-specific nucleotide gradients in euchromatin over evolutionary time. This orthogonality explains why GC-based zoning fails to predict mutation enrichment while Ev succeeds: the two statistics capture independent axes of sequence variation.

A natural concern is that Zone 3 mutation enrichment merely reflects transposable element density, since Alu SINEs are enriched 1.80-fold in Zone 3 relative to Zone 1 (RepeatMasker intersection, chromosome 1). However, the per-Alu-kilobase mutation rate is lower in Zone 3 than in Zone 1, indicating that Alu density and somatic mutation accumulation are independent co-consequences of heterochromatin state rather than causally linked. Both Alu insertion over evolutionary time and somatic C>T accumulation over a patient's lifetime reflect the same underlying property---constitutive heterochromatin with impaired surveillance---which Ev detects through the repeat-tract composition signature these elements leave in the sequence.

The replication timing analysis deserves particular attention. Late replication is often cited as a major confound in chromatin-mutation studies [4, 16], and it could in principle explain Zone 3 enrichment if Zone 3 windows were preferentially late-replicating. The stratified analysis (Table 4) rules this out: Zone 3 enrichment persists in all three replication timing strata, including in early-replicating regions where the Ev effect is actually strongest (OR = 1.37). This indicates that the sequence-composition signal captures chromatin-associated mutation rate variation that is not reducible to replication timing.

Several limitations should be noted. First, the formula parameters were generated with a fixed random seed (42) and were never optimized for any biological task. Of 20 random seeds tested, 12 produced qualitatively similar results (gene-dense regions scoring higher than pericentromeric regions), indicating moderate but not complete seed stability. Second, the ChIP-seq validation for H3K4me3 was limited to chromosome 1 due to data availability, and used a lymphoblastoid cell line (GM12878) rather than a breast-derived line. Cell-type matched data (e.g., MCF-7) would strengthen the mechanistic claims. Third, chromosome 19 represents a genuine edge case where the extreme gene density causes most of the chromosome to classify as Zone 3 by sequence composition despite being functionally euchromatic, producing an inverted OR. Excluding chromosome 19, all remaining 22 chromosomes show OR > 1.

The universal pan-cancer result (15/15 types, zero exceptions) and the constitutional germline baseline (Z3/Z1 = 1.13 in healthy individuals) together argue that Zone 3 vulnerability is a physical property of the genome encoded in its sequence composition. Cancer amplifies this baseline vulnerability 1.38-fold, likely through the global heterochromatin hypomethylation that characterizes tumor epigenomes [18]. The separation of tumor suppressors (Zone 3) from oncogenes (Zone 1) at OR = 9.2 provides a conceptual framework for understanding why certain driver genes are preferentially affected by accumulated mutations, and suggests that sequence composition alone can inform which genes are structurally predisposed to loss-of-function mutagenesis.

Looking ahead, this sequence-only approach may find applications in species where experimental epigenomic data are unavailable, in population-scale analyses where generating ChIP-seq for every sample is impractical, and in comparative genomics where one wishes to predict mutation-prone regions in newly assembled genomes. The finding that C>T transitions specifically drive Zone 3 enrichment also suggests that the approach might be adapted for predicting C>T mutation rates in clinical sequencing panels, where only targeted regions are interrogated.

The convergence of the cancer mutation signal with the epigenetic aging clock deserves particular comment. Horvath proposed that DNA methylation age reflects the cumulative work of an epigenetic maintenance system (EMS) [20], and observed that cancer tissues with high age acceleration paradoxically carry fewer somatic mutations, suggesting that the EMS protects against both aging-associated methylation drift and mutation accumulation. The present finding that aging clock CpGs concentrate in Zone 3 (OR = 1.94 for age-gaining CpGs) provides a mechanistic link: the same heterochromatic sequence environment that renders Zone 3 windows vulnerable to somatic mutation accumulation also renders the CpGs within them vulnerable to age-associated methylation gain when Polycomb-mediated maintenance declines. Cancer mutation geography and epigenetic aging thus appear to share a common origin in the sequence-encoded properties of heterochromatin. The fact that both can be predicted from DNA sequence alone, without experimental epigenomic data, suggests that the underlying chromatin vulnerability is a physical consequence of local sequence composition rather than a stochastic or purely regulatory phenomenon.

## Methods

### Ev formula

For each non-overlapping 5-kb window of genomic sequence, I counted the frequency of all 256 possible 4-mers. Windows with >10% ambiguous bases (N) were excluded. The 256-element frequency vector was normalized to sum to 1, then projected through a fixed random matrix P of dimensions 256 x 500, generated from a standard normal distribution using NumPy's default_rng(42) and divided by sqrt(256). From the 500-element projection, I selected the elements at the first 95 prime indices (indices 2, 3, 5, 7, ..., 499), computed the absolute skewness of these 95 values, and scaled the result: Ev = |skewness| x 6.07 + 0.10.

To remove the dominant GC-content effect, I computed the residual Ev_resid = (Ev - (-5.938 x GC + 4.471)) / 0.456, where the regression coefficients were derived from the genome-wide distribution of Ev values.

Zone thresholds were set at Ev_resid >= +0.382 (Zone 1) and Ev_resid <= -1.471 (Zone 3). These thresholds are fixed and were never optimized against any mutation data.

### Genome data

Human genome sequence (GRCh38/hg38) was obtained from UCSC Genome Browser. Chromosomes 1--22 and X were scanned in non-overlapping 5-kb windows, producing 582,028 classifiable windows. Zone assignments were cached as JSON files per chromosome.

### Somatic mutation data

Open-access masked somatic mutation MAF files were downloaded from the Genomic Data Commons (GDC) API (https://api.gdc.cancer.gov). For the primary breast cancer analysis, 992 TCGA-BRCA MAF files were obtained, yielding 89,565 genome-wide somatic mutations. For the pan-cancer analysis, 50 MAF files were randomly sampled from each of 15 TCGA projects (BRCA, LUAD, LUSC, COAD, PRAD, STAD, BLCA, LIHC, KIRC, HNSC, UCEC, THCA, SKCM, OV, GBM), totaling 750 patients. Each mutation was assigned to a zone based on the window it fell within. As an independent cross-validation, the MC3 pan-cancer mutation callset (3.6 million mutations across all TCGA projects) was used to verify the primary GDC-based results; all core findings replicated with the same direction and comparable effect sizes (MC3 BRCA OR = 1.585, 15/15 cancer types OR > 1).

### Statistical tests

For each analysis, I counted the number of mutations falling in Zone 1 and Zone 3 windows and computed a 2 x 2 contingency table (mutations present/absent x zone). Odds ratios and P-values were computed using Fisher's exact test (one-sided, testing Zone 3 > Zone 1). Mutation rates per window were computed as total mutations divided by total windows in each zone.

### ChIP-seq data

H3K9me3 narrowPeak data for GM12878 (all chromosomes) and H3K4me3 narrowPeak data for GM12878 (chromosome 1 only) were obtained from ENCODE. A 5-kb window was classified as "mark-positive" if any ChIP-seq peak overlapped it. Within each zone, I compared mutation rates in mark-positive versus mark-negative windows.

### Replication timing

GM12878 Repli-seq wavelet-smoothed signal (ENCFF001GVQ, bigWig format) was obtained from ENCODE. This dataset is aligned to hg19; at 5-kb window resolution, the coordinate offset from hg38 affects fewer than 0.5% of chromosome 1 windows. For each Zone 1 and Zone 3 window on chromosome 1, I computed the mean Repli-seq signal and divided windows into tertiles (early, mid, late). Zone-specific mutation rates and odds ratios were computed within each tertile.

### Germline data

High-coverage variant calls from the 1000 Genomes Project (2,504 individuals) were obtained for chromosomes 1--5 from the EBI FTP server. The 26.2 million germline variants were intersected with the zone map using the same 5-kb window framework.

### Methylation data

Illumina 450K methylation array data for 32 matched tumor-normal pairs from TCGA-BRCA were obtained from GDC. For each pair, mean methylation (beta-values) was computed for probes falling within Zone 1 and Zone 3 windows. The difference in methylation between tumor and normal (Delta-beta) was tested by paired t-test.

### Cross-species analysis

Chromosome-scale FASTA sequences for ten eukaryotic species were obtained from NCBI RefSeq and Ensembl (see Additional file 1 for accession numbers). Each genome was scanned with the identical Ev formula and zone thresholds.

### Drug target classification

101 Tier 1 genes from the COSMIC Cancer Gene Census were classified by computing the Ev zone for each 5-kb window overlapping the gene body. A gene was assigned to its majority zone (the zone containing >50% of its windows). The proportion of tumor suppressors versus oncogenes in each zone was compared using Fisher's exact test (two-sided).

### Epigenetic aging clock analysis

The 353 CpG sites composing Horvath's multi-tissue age predictor were obtained from the supplementary data of [20] (Additional file 3 of the original publication). Genomic coordinates (hg18/Build 36) were mapped to the hg38 zone framework by matching each CpG to the 5-kb window containing its position, with a tolerance of +/-2 windows (+/-10 kb) to accommodate coordinate system differences; at 5-kb resolution, this tolerance captures >95% of CpGs correctly. Each mapped CpG was classified into Zone 1, 2, or 3 based on the Ev_resid of its containing window. CpGs were further stratified by their marginal age relationship (positive = hypermethylated with age, n = 193; negative = hypomethylated with age, n = 160) as reported in the original publication. Zone enrichment was tested against the genome-wide zone distribution (582,028 windows) using Fisher's exact test and chi-squared tests.

### Reverse-mapping and repeat element analysis

Pearson partial correlations between individual 4-mer frequencies and Ev_resid, controlling for GC content, were computed across all valid chromosome 1 windows (n = 46,105). The L2 norm of each 4-mer's weight vector through the projection matrix was computed to test for systematic projection bias between CpG-containing and non-CpG 4-mers. RepeatMasker annotations (hg38) were downloaded from UCSC Table Browser. For each 5-kb window on chromosome 1, total base pairs covered by each transposable element family were computed and correlated with Ev_resid and mutation counts. Per-Alu mutation rates were computed as total mutations divided by total Alu kilobases within each zone. Principal component analysis of the 256-dimensional 4-mer frequency matrix was performed using singular value decomposition; AT-skew was defined as (A-T)/(A+T) computed per window. Partial correlations with Ev_resid after GC removal were computed using ordinary least squares residualization.

### Software and reproducibility

All analyses were performed in Python 3 using NumPy, SciPy, and pyBigWig. The Ev formula uses a deterministic random matrix generated by numpy.random.default_rng(42); the identity of the matrix was verified by asserting P[0,0] = 0.01904482. Code is available at [repository URL upon publication].

## Declarations

### Ethics approval and consent to participate
Not applicable. All data used in this study are publicly available through open-access repositories (GDC, ENCODE, 1000 Genomes Project).

### Availability of data and materials
All data sources are listed in Additional file 1. The Ev formula implementation and analysis scripts will be deposited in a public repository upon publication.

### Competing interests
The author declares no competing interests.

### Funding
This research was conducted independently without external funding.

### Author contributions
Aditya Tiwari conceived the study, developed the Ev formula, performed all analyses, and wrote the manuscript.

### Acknowledgement
AI assistance (Claude, Anthropic) was used for code development and manuscript editing. All scientific analyses, experimental design, data interpretation, and conclusions are the author's own.

## References

1. Schuster-Bockler B, Lehner B. Chromatin organization is a major influence on regional mutation rates in human cancer cells. Nature. 2012;488:504--507.

2. Pleasance ED, Cheetham RK, Stephens PJ, et al. A comprehensive catalogue of somatic mutations from a human cancer genome. Nature. 2010;463:191--196.

3. Makova KD, Hardison RC. The effects of chromatin organization on variation in mutation rates in the genome. Nat Rev Genet. 2015;16:213--223.

4. Stamatoyannopoulos JA, Adzhubei I, Thurman RE, et al. Human mutation rate associated with DNA replication timing. Nat Genet. 2009;41:393--395.

5. Polak P, Karlic R, Koren A, et al. Cell-of-origin chromatin organization shapes the mutational landscape of cancer. Nature. 2015;518:360--364.

6. Tolstorukov MY, Volfovsky N, Stephens RM, Park PJ. Impact of chromatin structure on sequence variability in the human genome. Nat Struct Mol Biol. 2011;18:510--515.

7. Akdemir KC, Le VT, Chandran S, et al. Disruption of chromatin folding domains by somatic genomic rearrangements in human cancer. Nat Genet. 2020;52:294--305.

8. Supek F, Lehner B. Differential DNA mismatch repair underlies mutation rate variation across the human genome. Nature. 2015;521:81--84.

9. Liu C, Wang Z, Wang J, et al. Predicting regional somatic mutation rates using DNA motifs. PLoS Comput Biol. 2023;19:e1011536.

10. Alexandrov LB, Kim J, Haradhvala NJ, et al. The repertoire of mutational signatures in human cancer. Nature. 2020;578:94--101.

11. Tomkova M, Tomek J, Kriaucionis S, Schuster-Bockler B. Mutational signature distribution varies with DNA replication timing and strand asymmetry. Genome Biol. 2018;19:129.

12. Fryxell KJ, Moon W-J. CpG mutation rates in the human genome are highly dependent on local GC content. Mol Biol Evol. 2005;22:650--658.

13. Supek F, Lehner B. Clustered mutation signatures reveal that error-prone DNA repair targets mutations to active genes. Cell. 2017;170:534--547.

14. Haradhvala NJ, Polak P, Stojanov P, et al. Mutational strand asymmetries in cancer genomes reveal mechanisms of DNA damage and repair. Cell. 2016;164:538--549.

15. Jinks-Robertson S, Bhagwat AS. Transcription-associated mutagenesis. Annu Rev Genet. 2014;48:341--359.

16. Koren A, Polak P, Nemesh J, et al. Differential relationship of DNA replication timing to different forms of human mutation and variation. Am J Hum Genet. 2012;91:1033--1040.

17. Berman BP, Weisenberger DJ, Aman JF, et al. Regions of focal DNA hypermethylation and long-range hypomethylation in colorectal cancer coincide with nuclear lamina-associated domains. Nat Genet. 2012;44:40--46.

18. Timp W, Feinberg AP. Cancer as a dysregulated epigenome allowing cellular growth advantage at the expense of the host. Nat Rev Cancer. 2013;13:497--510.

19. Bernardi G. Isochores and the evolutionary genomics of vertebrates. Gene. 2000;241:3--17.

20. Horvath S. DNA methylation age of human tissues and cell types. Genome Biol. 2013;14:R115.

## Figure legends

**Figure 1. Ev formula and zone classification.** (a) Schematic of the Ev computation: raw DNA sequence to 4-mer frequencies to random projection to skewness to zone classification. (b) Genome-wide distribution of Ev_resid across 582,028 windows, with Zone 1 (blue, Ev_resid >= +0.382) and Zone 3 (red, Ev_resid <= -1.471) thresholds marked. (c) Karyotype-style map of Zone 1 (blue) and Zone 3 (red) assignments across human chromosomes 1--22 and X.

**Figure 2. Per-chromosome odds ratios.** Forest plot showing Zone 3 vs Zone 1 mutation enrichment odds ratios and 95% confidence intervals for each chromosome. Dashed vertical line at OR = 1.0; genome-wide OR = 1.682 (gray band). Chromosome 19 is the sole inversion (OR = 0.90); chromosome X shows the strongest effect (OR = 2.24).

**Figure 3. Pan-cancer universality.** Forest plot of Zone 3 mutation enrichment odds ratios across 15 TCGA cancer types. All 15 show OR > 1 with no exceptions. PRAD (prostate) and UCEC (uterine) show the highest enrichment; GBM (brain) and LUSC (lung squamous) the lowest.

**Figure 4. C>T specificity.** Bar chart showing Zone 3 / Zone 1 odds ratios for each of six substitution classes in TCGA-BRCA. C>T is the only class with OR > 1 (OR = 1.24, P = 8.7 x 10^-20^). T-origin mutations are depleted in Zone 3.

**Figure 5. Replication timing independence.** (a) Zone 3 vs Zone 1 odds ratios within each replication timing tertile (early, mid, late). Zone 3 enrichment is significant in all three strata. (b) Cross-tabulation showing approximately equal numbers of Zone 3 windows in each replication timing class, confirming that Zone 3 is not a proxy for late replication.

**Figure 6. Constitutional vulnerability and cancer amplification.** (a) Germline variant density (1000 Genomes, 2,504 individuals) across zones: Z1 < Z2 < Z3 (P = 1.5 x 10^-82^). (b) Comparison of germline Z3/Z1 ratio (1.13) with somatic Z3/Z1 ratio (1.55), showing 1.38-fold cancer amplification. (c) Methylation change (Delta-beta) in 32 matched tumor-normal pairs: Zone 3 hypomethylation (Delta-beta = -0.054) and Zone 1 hypermethylation (Delta-beta = +0.032).

**Figure 7. Epigenetic aging clock CpGs concentrate in Zone 3.** (a) Distribution of 353 Horvath clock CpGs across Ev zones compared to genome-wide expectation. Clock CpGs are depleted from Zone 1 (OR = 0.53, P = 5.5 x 10^-6^) and enriched in Zone 3 (OR = 1.37, P = 8.4 x 10^-3^). (b) Stratification by age-relationship direction: the 193 positive (age-gaining) clock CpGs show 1.94-fold Zone 3 enrichment (P = 2.0 x 10^-5^), while the 160 negative (age-losing) CpGs concentrate in Zone 2 with no Zone 3 enrichment. (c) Schematic illustrating the convergence of cancer mutation enrichment and epigenetic aging in Zone 3 heterochromatin.

## Additional files

**Additional file 1.** Supplementary tables and figures, including: complete per-chromosome results, pan-cancer signature analysis (LUAD, SKCM, PRAD), drug target gene classification with TSG/oncogene segregation statistics, ChIP-seq enrichment details, cross-species accession numbers, seed stability analysis (12/20 seeds), TMB stratification, reverse-mapping and repeat element analysis, MC3 cross-validation results, and survival analysis (negative result).

## Tables

**Table 1. Genome-wide Zone 3 mutation enrichment in TCGA-BRCA.**

| | Zone 1 | Zone 2 | Zone 3 |
|---|---|---|---|
| Windows | 152,519 | 287,733 | 141,776 |
| Mutations | 18,413 | 44,547 | 26,605 |
| Rate (mutations/window) | 0.121 | 0.155 | 0.188 |
| Z3/Z1 ratio | | | 1.554 |
| Fisher OR (Z3 vs Z1) | | | **1.682** |
| P-value | | | **< 10^-300^** |

992 patients, 89,565 total mutations, 582,028 windows.

**Table 2. Pan-cancer Zone 3 mutation enrichment across 15 TCGA cancer types.**

| Cancer type | n patients | OR | P |
|---|---|---|---|
| TCGA-BRCA (breast) | 50 | 1.53 | 9.5 x 10^-14^ |
| TCGA-LUAD (lung adeno) | 50 | 1.52 | 1.6 x 10^-70^ |
| TCGA-LUSC (lung squam) | 50 | 1.43 | 6.2 x 10^-54^ |
| TCGA-COAD (colon) | 50 | 1.57 | 2.3 x 10^-157^ |
| TCGA-PRAD (prostate) | 50 | 1.91 | 4.0 x 10^-105^ |
| TCGA-STAD (stomach) | 50 | 1.65 | 1.1 x 10^-145^ |
| TCGA-BLCA (bladder) | 50 | 1.52 | 5.0 x 10^-74^ |
| TCGA-LIHC (liver) | 50 | 1.67 | 1.5 x 10^-52^ |
| TCGA-KIRC (kidney) | 50 | 1.65 | 3.1 x 10^-21^ |
| TCGA-HNSC (head/neck) | 50 | 1.69 | 3.4 x 10^-99^ |
| TCGA-UCEC (uterine) | 50 | 1.80 | < 10^-300^ |
| TCGA-THCA (thyroid) | 50 | 1.57 | 4.5 x 10^-5^ |
| TCGA-SKCM (melanoma) | 50 | 1.77 | 1.6 x 10^-204^ |
| TCGA-OV (ovarian) | 50 | 1.74 | 4.4 x 10^-38^ |
| TCGA-GBM (brain) | 50 | 1.45 | 2.0 x 10^-44^ |

Mean OR = 1.63. Zero exceptions to Z3 > Z1.

**Table 3. Mutation class-specific Zone 3 enrichment in TCGA-BRCA.**

| Substitution class | Z3/Z1 ratio | OR | P |
|---|---|---|---|
| **C>T** | **1.120** | **1.243** | **8.7 x 10^-20^** |
| C>A | 0.972 | 0.968 | 0.33 |
| C>G | 0.975 | 0.969 | 0.30 |
| T>A | 0.784 | 0.775 | 7.9 x 10^-6^ |
| T>C | 0.760 | 0.742 | 1.8 x 10^-11^ |
| T>G | 0.713 | 0.703 | 7.1 x 10^-9^ |

**Table 4. Zone 3 enrichment stratified by replication timing (chromosome 1).**

| RT stratum | Z1 windows | Z3 windows | OR | P |
|---|---|---|---|---|
| Early | 3,213 | 3,327 | 1.372 | 5.5 x 10^-5^ |
| Mid | 3,078 | 3,348 | 1.284 | 6.8 x 10^-4^ |
| Late | 3,248 | 3,245 | 1.157 | 3.8 x 10^-2^ |

39,059 windows with both zone assignment and replication timing data. GM12878 Repli-seq (ENCFF001GVQ).
