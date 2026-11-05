# Methods — Proteomics cercevesine gore yeniden yazim
Tarih: 2026-11-04 | A1 ve A2 taslaklarinin devami

## BIRINCIL ANALIZ KARARI

Ters-varyans agirliklama BIRAKILDI. Birincil analiz agirliksiz; ters-varyans
ve kirpilmis semalar duyarlilik olarak raporlanacak. Gerekce Methods 2.3'te.

Agirliksiz sonuc: 4/4 kosulda permutasyon kontrolu geciliyor
(P = 0.030, 0.005, 0.020, 0.005). Ters-varyans ile 2/4.

---

# TASLAK METIN

## 2.1 Scope and reuse design

The study is a secondary computational analysis of publicly deposited
proteomics data. No new mass-spectrometry, biochemical, cell, animal or human
experiments were performed, and no raw spectra were reprocessed: the deposited
dataset provides processed site-level quantification, and we used it as
deposited. The claim boundary was fixed before analysis. We asked whether local
structural descriptors carry reproducible information about site-specific CML
reduction within bovine serum albumin, and we treated condition-level
replication as the primary test. We make no claim about proteome-wide
repairability, about enzyme delivery, or about transfer to other protein
scaffolds.

## 2.2 Data source and condition identification

Site-level quantification was taken from ProteomeXchange dataset PXD070124,
files `BSA_{0p25,0p5,1p0,2p0}_865_794_processed_filt_quant.csv`. SHA-256
prefixes of the four files are reported in the data availability statement. Each
file reports, per lysine, the protein, residue number, modification, peptide
sequence, an identification confidence score, and percent abundance for three
active-enzyme and three inactive-control replicates. The four files correspond
to the four substrate loads described in the source report, in which 250, 500,
1000 or 2000 ng of CML-modified BSA in a 25 uL reaction was incubated with
enzyme at fixed concentration; the enzyme-to-substrate ratio therefore falls as
load rises. Residue numbering in the deposited files is UniProt numbering for
P02769, which we confirmed against the mature-chain offset of 24 used in the
source analysis.

The source report does not state which substrate load its site-level figure
represents. We identified it as the 500 ng condition on three independent
grounds, reported in Results: exact agreement of replicate-level percent
abundance values with the figure source data, agreement of the number of sites
passing the confidence threshold, and agreement of the number of sites showing
more than ninety percent reduction. Sites were retained when the identification
confidence score exceeded 90, the threshold used in the source analysis.

## 2.3 Site-level response labels and the weighting decision

For each site in each condition the response was defined as one minus the ratio
of the active-enzyme mean percent abundance to the inactive-control mean percent
abundance, computed from the three replicates on each arm and truncated to the
unit interval. Sites with a control mean of zero were excluded because the ratio
is undefined. Truncation affected a small number of sites per condition (two,
one, five and two at 250, 500, 1000 and 2000 ng), and the number of sites at the
lower and upper boundary after truncation is reported in Results. The endpoint is
a relative reduction in CML signal under one assay and is not interpreted as an
intrinsic probability that a site is repairable.

An approximate variance was obtained for each response by the delta method,
combining the replicate-mean variances of the two arms as
var(A)/C^2 + A^2 var(C)/C^4, where A and C are the active and control means. The
source pipeline converted these variances into inverse-variance sample weights.
We evaluated that choice rather than inheriting it, and we did not retain it.

Three observations motivate the decision. First, the weights are extremely
skewed: within a condition they span roughly three to five orders of magnitude,
and the Kish effective sample size falls to between 0.23 and 0.32 of the nominal
site count, so the thirty-three nominal sites of the published condition
correspond to an effective sample of 9.7. Second, the weighting suppresses the
structure-response relationship: under inverse-variance weighting the
permutation-controlled signal was recovered in three of four conditions, whereas
without weighting it was recovered in all four, with higher rank correlations in
three of the four. Third, clipping the weights to a moderate range reproduced the
weighted outcome rather than the unweighted one, which indicates that the problem is
the direction of the weighting rather than a few extreme values.

The mechanism is specific to ratio-based response labels. The delta-method
variance carries the control mean in the denominator, so sites with low control
abundance receive large variances and correspondingly small weights. Low
occupancy, however, is a property of the site rather than noise about the
structure-response relationship, so inverse-variance weighting systematically
downweights low-occupancy sites and introduces an occupancy-dependent bias. We
therefore report unweighted analyses as primary and inverse-variance and clipped
weighting as sensitivity analyses, and we recommend that reuse of ratio-based PTM
response labels treat inverse-variance weighting as a modelling choice requiring
justification rather than a default.

## 2.4 Structural descriptors

Structural features were extracted from the crystallographic model of bovine
serum albumin, PDB 4F5S, chain A. Solvent-accessible surface area was computed
with the Shrake-Rupley algorithm and normalised to the residue-specific maximum
accessibility scale for lysine. The representation contains nine descriptors:
residue solvent-accessible surface area, lysine-relative relative accessibility,
distance to the protein centroid, hydrophobic-residue fraction, cysteine count
and aromatic-residue count within 8 A, and C-alpha contact counts within 8, 10
and 12 A. Contact counts are treated as static packing descriptors and not as
measures of flexibility.

Two feature classes were excluded deliberately. Crystallographic B-factors,
available for this scaffold, were excluded because they have no equivalent in
predicted structures and their inclusion would make the representation
scaffold-specific. Three contact-derived variables previously proposed as
flexibility proxies were excluded because their interpretation as dynamics is
unvalidated and they did not improve performance in ablation. C-alpha
coordinates were used to define spatial validation groups and were never
supplied as predictive features.

Descriptors were extracted for every lysine appearing in the union across the
four conditions. We note that the descriptor table distributed with the original
analysis covers only the thirty-three sites of the published condition, so eleven
lysines appearing at other substrate loads (positions 36, 44, 151, 160, 228, 235,
248, 285, 340, 544 and 548) were absent from it. Their absence reflects the
single-condition scope of that table rather than any structural or quality
criterion: all eleven are ordinary lysines in chain A of 4F5S and yield
descriptors under the rules above. We therefore extracted the full set and
verified the extraction by recomputing descriptors for the thirty-three
previously covered sites, which reproduced the distributed values for all
thirteen fields to floating-point precision (maximum absolute difference
1.4e-14).

## 2.5 Spatially blocked validation and model

Because neighbouring lysines share local environment, random site-level
cross-validation would place chemically similar residues in both training and
test folds. Sites were therefore partitioned by their C-alpha coordinates using
Ward agglomerative clustering, into four groups for full condition sets and three
groups for the twenty-site analyses, and all members of a group were withheld
together. Standardisation was fitted on training folds only.

The regressor was an extremely randomised trees ensemble with 500 estimators,
maximum depth 4, minimum leaf size 2 and a fixed random seed, applied identically
in every condition and every analysis. The null comparator was the training-fold
mean response. Performance is summarised by mean absolute error and Spearman rank
correlation, with relative improvement over the mean predictor. We do not report
the coefficient of determination, which is unreliable for model assessment at
this sample size.

## 2.6 Site-structure permutation control

The primary inferential test asks whether prediction depends on each response
being paired with its own local environment. Complete structural feature vectors
were permuted among sites while labels, spatial partitions and the feature
covariance structure were held fixed, and the full blocked analysis was repeated
on each permutation. Two hundred permutations were used per analysis, and the
empirical P value is the proportion of permutations whose mean absolute error was
no greater than the observed value, with the usual plus-one correction. Because
the permutation preserves the multivariate feature distribution and the model
architecture, a significant result cannot be explained by generic properties of
the feature matrix.

## 2.7 Matched-site core and random-subset reference

Twenty lysines passed the confidence threshold in all four conditions and also
had extractable descriptors. Restricting every condition to this matched set
removes site composition as a difference between conditions, at the cost of
reducing the sample. To separate loss of power from site composition we drew two
hundred random twenty-site subsets of the published condition without
replacement, repeated the blocked analysis on each, and located the matched set
within the resulting distribution of rank correlations and baseline improvements.
We note that the random subsets are drawn from the same thirty-three sites and
therefore overlap the matched set substantially, so the reference distribution is
not independent of it; the comparison bounds the power explanation but does not
establish that the matched sites form a distinct class.

## 2.8 Confound tests

To test whether structural descriptors predict measurability rather than
response, we related each descriptor to the number of conditions in which a site
passed the confidence threshold and to its mean confidence score, using Spearman
correlations with a Bonferroni threshold of 0.05/9. We compared matched and
unmatched site groups descriptor by descriptor with Mann-Whitney tests under the
same threshold, tested the association between confidence score and response
within the published condition, and compared the response distributions of the
two site groups.

## 2.9 Protein-level comparison

The frozen model was compared with the bulk CML-ELISA responses reported for
casein, haemoglobin and collagen in the source study, using canonical bovine
components with stoichiometric weighting because the commercial species and
isoforms were not specified. The observed side was read at full resolution: the
source report deposits three modification levels per scaffold, and the ordering
was evaluated within each level rather than after collapsing levels to a scaffold
median, since a three-point rank correlation cannot attain a P value below one in
six. Predicted scaffold scores are level-invariant, so the independent units for
the model-observation comparison are the three scaffolds and the attainable P
value floor for an exact-ordering match is one in six regardless of how many
levels are measured. The comparison was prespecified as an applicability stress
test rather than site-resolved external validation, because the endpoint is a
bulk immunoassay that does not identify which lysines carried CML.

## 2.10 Software and reproducibility

Analyses were implemented in Python using scikit-learn, SciPy, NumPy, pandas and
Biopython, the last for parsing coordinates and computing solvent accessibility.
All random seeds, feature names, partition assignments and input and output
checksums were recorded. Code, processed label tables, configurations and result
files are available at the repository cited in the data availability statement.
The four source files are identified there by SHA-256 so that the derived label
tables can be regenerated from the archive without reprocessing raw spectra.

---

## DEGISMESI GEREKEN YERLER (A1 ve A2 taslaklarinda)

Agirliklama kararindan sonra su bolumler revize olacak:

### A1 Tartisma
- "matched-site analysis is harder to explain away" paragrafi: artik agirliksiz
  sonuclarla yazilacak (cekirdek yine basarisiz, percentile %3)
- 2/4 tekrar anlatisi 4/4 olarak degisecek — Tartisma'nin ikinci paragrafi
  tamamen yeniden yazilacak
- Agirliklama bulgusu YENI bir paragraf olarak girecek (metodolojik katki)

### A1 Ozet
- "two of four conditions" -> "all four conditions (unweighted)"
- agirliklama bulgusu bir cumle olarak eklenecek
- ESS bilgisi eklenecek

### A2 Results
- 3.3 yeniden yazilacak: 4/4 agirliksiz, 2/4 ters-varyans
- 3.4 (dinamik aralik) korunuyor ama artik farkli bir soruya cevap veriyor
- 3.5 agirliksiz sayilarla guncellenecek (cekirdek P = 0.080/0.194/0.264/0.129,
  percentile %3)
- YENI 3.x: agirliklama duyarliligi ve ESS
- Tablo 2'ye agirliksiz kolonlar eklenecek
- Sekil 2'ye agirliklama karsilastirma paneli eklenecek
