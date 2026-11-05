# JPR icin revize manuscript — Title / Abstract / Introduction / Results / Discussion
Tarih: 2026-11-04 | Nihai 44-bolge sayilariyla

## JPR CERCEVESININ MANTIGI

JPR kapsam maddesi: "Computational methods including bioinformatics, artificial
intelligence and machine learning, or other software tools to aid in data
analysis and analysis of experimental systems."

Bu yuzden vurgu sirasi DEGISTI. Proteomics taslaginda yapi-yanit sinyali one
cikiyordu; JPR'de **etiket kurulumu ve veri yeniden kullanimi** one cikiyor,
yapi-yanit analizi ise bunun gosterildigi arac oluyor. Bu hem kapsama daha
uygun hem de elimizdeki en guclu bulgu artik gercekten o.

---

# TASLAK METIN

## Baslik

**Condition identity and label weighting govern structure-response inference
from reused site-resolved PTM data: a reanalysis of enzymatic
carboxymethyllysine repair**

## Ozet

Site-resolved measurements of how a post-translational modification responds to
an enzyme are increasingly reused as labels for computational models, but little
is known about how robust such labels are to the analytical choices made when
deriving them. We examined this using ProteomeXchange PXD070124, which
accompanies the first report of an enzyme that reverses
N-epsilon-carboxymethyllysine and deposits site-level quantification for four
substrate loads. The published site-level figure reports one load without stating
which; we identify it as the 500 ng condition from replicate-level agreement, the
number of sites passing the confidence threshold, and the number of highly
responsive sites. Converting active-versus-control abundance ratios into
inverse-variance sample weights, a natural choice given the available
delta-method variance, proved harmful: the weights spanned up to five orders of
magnitude, the Kish effective sample size fell to 0.23-0.32 of the nominal site
count, and a structure-response relationship was recovered in three of four
conditions under weighting against four of four without it. Clipping the weights
did not repair this. The mechanism is that the delta-method variance carries the
control mean in its denominator, so low-occupancy sites are systematically
downweighted. Unweighted, local structural descriptors carried
permutation-controlled predictive information in all four loads (P = 0.010,
0.005, 0.005, 0.005). The relationship was not uniform across sites: the twenty
lysines quantified confidently at every load showed no signal in any condition,
at the third percentile of random-subset resampling, which excludes loss of
power. Reuse of ratio-based site-level response data therefore requires that the
reporting condition be stated and that inverse-variance weighting be justified
rather than assumed.

## Giris

Site-resolved proteomics increasingly supplies not only maps of where a
post-translational modification occurs, but quantitative measurements of how
individual modified residues respond to a perturbation. Such measurements are
attractive as labels for computational models, because they promise to turn a
chemical modification into a per-residue phenotype. They are also scarce, small,
and produced under conditions chosen for the original experiment rather than for
reuse. How robust they are to the analytical choices a secondary user must make
has not, to our knowledge, been examined on a real dataset.

Prediction of lysine glycation from sequence and structure is itself an
established line of work, running from early neural-network classifiers through
structure-aware support vector machines to recent transformer-based models.
Every one of these predicts occurrence: which lysines carry the modification.
None predicts how a modified lysine responds to an enzyme, because until
recently no such measurements existed. Response labels are therefore a new and
largely unexamined class, and the question of how robustly they can be derived
from a deposition has not been asked.

N-epsilon-carboxymethyllysine (CML) provides an unusually clean test case.
Chemically stable and abundant in long-lived proteins, CML was treated as a
terminal modification: endogenous defences such as the glyoxalase system and
fructosamine-3-kinase act on precursors or early glycation adducts, and the
reported deglycase activity of DJ-1/PARK7 remains contested. The recent
description of an engineered FAD-dependent oxidase that oxidises the
N-epsilon-carboxymethyl group and restores native lysine therefore created the
first site-level response measurements for this modification class, together with
a public deposition.

The source study reported that of thirty-three CML-modified lysines detected on
bovine serum albumin, thirty showed reduced CML, and it noted qualitative
structural trends: activity was "not strictly a function of surface
accessibility", three lysines in surface-exposed helices "showed varied CML
reductions of 90%, 38%, and 0%", activity "weakly correlated with local protein
flexibility (B-factor)" such that "flexibility may be necessary but is not
sufficient for catalysis", and "some local sequence preference was observed, with
CMLase showing higher activity on sites with neighboring cysteines" while
"reduced enzymatic activity was associated with bulky hydrophobic or aromatic
residues". The legend of the corresponding figure nonetheless states that the
enzyme "has varying activity across all sites with no clear preference based on
primary or secondary structure". These observations are the starting point for
the present work rather than a gap in it: the structural hypotheses have been
stated, and what is missing is a quantitative test of whether they survive the
choices a reuser must make.

Those choices turn out to be consequential, and they are the subject of this
paper. The deposited dataset contains processed site-level quantification for
four substrate loads, whereas the published site-level figure reports one of them
without identifying it. A response defined as one minus the ratio of
active-enzyme to inactive-control abundance has a readily computed delta-method
variance, which invites inverse-variance weighting. And the set of lysines that
pass an identification-confidence threshold differs between loads, so the site
composition of a label set is itself a function of the condition from which it
was taken.

We therefore treat the reuse operation as the object of study. We identify which
condition the published figure reports, reproduce it with an independent
implementation, and then ask how the apparent structure-response relationship
depends on the weighting scheme, on the substrate load, and on which sites are
included. The structural model is a vehicle for these questions rather than a
contribution in itself; the descriptors and the regressor are standard, and we
make no claim about proteome-wide repairability or transfer to other scaffolds.

## Sonuclar (Results)

### Condition identity

The deposited dataset provides processed site-level quantification for four
substrate loads, each with three active-enzyme and three inactive-control
replicates (Table 1); the overall workflow is summarised in Figure 1. The source report does not state which load its site-level
figure represents. Three independent lines of evidence identify it as the 500 ng
condition. The replicate-level percent-abundance values underlying the published
figure match that condition exactly; for example the lysine at position 65
carries active replicates of 9.98, 8.63 and 9.43 and inactive replicates of
17.64, 17.40 and 17.92 in both. The 500 ng condition is the only one in which
exactly thirty-three lysines pass the confidence threshold used in the source
analysis, the other loads yielding twenty-two, thirty-seven and thirty-eight,
and the source report describes thirty-three detected sites. Seven sites in that
condition show more than ninety percent reduction, the number reported in the
source study.

Two counts differ slightly and we record rather than reconcile them. We find
thirty-two of thirty-three sites with any reduction where the source reports
thirty, and twenty sites above fifty percent reduction where the source reports
twenty-one. Both are consistent with different handling of values at the
truncation boundary and neither affects condition identification, which rests on
the exact replicate match.

### An independent implementation reproduces the published condition

Applied to the 500 ng condition under the weighting scheme of the source
pipeline, the nine-descriptor representation reached a mean absolute error of
0.184 against a mean-predictor baseline of 0.300, with Spearman rank correlation
0.681 and a permutation-control advantage of 0.100 (P = 0.005). The corresponding
published values for the same sites are 0.186, 0.699 and 0.113. The agreement
establishes that the feature definitions, spatial blocking and permutation
control reproduce the published behaviour, so the results that follow are not
attributable to implementation differences.

### Label weighting determines whether the relationship is detected

Inverse-variance weights derived from the delta-method variance are extremely
skewed. Within a condition they span five orders of magnitude, from 0.0002 to
28.4 at the lowest load, and the Kish effective sample size falls to between 0.23
and 0.32 of the nominal site count: the thirty-three nominal sites of the
published condition correspond to an effective sample of 9.7 (Table 2,
Figure 2A and 2B).

The consequence is substantive. Under inverse-variance weighting the
permutation-controlled signal was recovered in three of four conditions, failing
at the lowest load (Spearman 0.141, P = 0.109). Without weighting it was
recovered in all four, and rank correlations were higher in three of the four
conditions (Table 3, Figure 2C and 2D). Clipping the weights to a moderate range
reproduced the weighted result rather than the unweighted one, which indicates that the
difficulty is the direction of the weighting rather than a few extreme values.

The mechanism is specific to ratio-based response measures. The delta-method
variance carries the control mean in its denominator, so sites with low control
abundance receive large variances and correspondingly small weights. Low
occupancy is a property of the site rather than noise about the
structure-response relationship, so inverse-variance weighting systematically
downweights low-occupancy sites. We therefore report unweighted analyses as
primary.

### The relationship replicates across all four substrate loads

Unweighted, local structural descriptors carried permutation-controlled
predictive information in every condition: empirical P values were 0.010, 0.005,
0.005 and 0.005 at 250, 500, 1000 and 2000 ng, with Spearman correlations of
0.414, 0.616, 0.399 and 0.581 and improvements over the mean predictor of 0.090,
0.102, 0.064 and 0.091 (Table 3). Three of the four reach the resolution limit of
two hundred permutations; the observed errors lie below their permutation nulls in
every condition (Figure 3A, 3B).

Extending the descriptor table to every lysine in the union across conditions
contributed to this result. The table distributed with the original analysis
covers only the thirty-three sites of the published condition; eleven lysines
appearing at other loads were absent from it, not because of any structural or
quality criterion but because the table was built for one condition. Including
them raised the 1000 ng condition from P = 0.020 to P = 0.005 and the 250 ng
condition from P = 0.030 to P = 0.010.

### Response dispersion does not explain the differences between conditions

A natural explanation for a relationship appearing more strongly in some
conditions than others is that response variation collapses in the weaker ones,
leaving less to predict. The data do not support it. Response standard
deviations were 0.357, 0.310, 0.304 and 0.298 across the four loads and
interquartile ranges 0.574, 0.509, 0.533 and 0.492, a spread of at most twenty
percent between the widest and narrowest condition. Mean response declined with
increasing substrate load, from 0.602 at 500 ng to 0.455 at 2000 ng, which is the
expected direction when enzyme concentration is held constant and substrate
rises, and serves as a sanity check that the conditions behave as the design
implies. The conditions differed mainly in where responses accumulated at the
upper boundary, with seven sites above 0.90 at 500 ng against five, zero and one
elsewhere, but total dispersion was comparable throughout.

### Site composition, not power, limits the relationship

Twenty lysines passed the confidence threshold in all four conditions.
Restricting every condition to this matched set removes site composition as a
difference between conditions, and under that restriction no condition retained a
permutation-controlled signal: P values were 0.080, 0.194, 0.264 and 0.129, and
in the published condition the rank correlation became negative at -0.224
(Table 4).

Because reducing the sample also reduces power, we drew two hundred random
twenty-site subsets of the published condition. Their mean Spearman correlation
was 0.400, the median 0.461, and seventy-four percent exceeded 0.3; mean
improvement over baseline was 0.062. Twenty sites therefore suffice to recover
the relationship when drawn at random. The matched set reached -0.224 and a
baseline improvement of 0.003, placing it at the third percentile (Figure 3C, 3D). Power is thus
excluded as the explanation. We note that the random subsets are drawn from the
same thirty-three sites and overlap the matched set substantially, so the
reference distribution is not independent of it, and a third-percentile outcome
arises by chance in roughly one draw in thirty; the comparison bounds the power
explanation without establishing that the matched sites form a distinct class.

### Measurement confidence does not confound the relationship

Structural descriptors did not predict in how many conditions a site was
confidently quantified after correction for testing nine descriptors (Bonferroni
threshold 0.0056); the strongest association was neighbourhood hydrophobic
fraction (Spearman 0.442, P = 0.0101), which does not survive correction, and no
descriptor was associated with mean confidence (all P > 0.11). Within the
published condition, confidence score and response were unrelated (Spearman
-0.087, P = 0.63). The matched and unmatched site groups had indistinguishable
response distributions, with means of 0.602 and 0.603 (Mann-Whitney P = 0.84).
Comparing the groups descriptor by descriptor gave no difference surviving
correction, the smallest P value again belonging to neighbourhood hydrophobic
fraction (0.404 against 0.251, P = 0.0164). We report the convergence of
hydrophobic neighbourhood across two uncorrected tests as a lead, noting that it
is the descriptor family the source study associated with reduced activity, and
not as a finding.

## Tartisma (Discussion)

The object of this study is not an enzyme but a class of data. Site-resolved
measurements of how a modification responds to a perturbation are beginning to
be reused as quantitative labels, and the steps between an archive and a
conclusion are numerous enough that the same deposition can support different
answers. Working through those steps on one dataset, we found two that mattered
enough to change the result, and both were encountered rather than anticipated.

The first and more general is label weighting. Defining a response as one minus
the ratio of treated to control abundance makes a delta-method variance
immediately available, and using it as an inverse-variance weight is the
textbook response to heteroscedastic measurements. Here it was
counterproductive. The weights spanned up to five orders of magnitude, the effective
sample size fell to roughly a third of the nominal site count, and the
structure-response relationship was recovered in three of four conditions
weighted against four of four unweighted, with lower rank correlations in three
of the four. That clipping the weights reproduced the weighted rather than the
unweighted outcome locates the problem in the direction of the weighting, not in
outliers. The mechanism is worth stating plainly because it is not specific to
this dataset: the delta-method variance of a ratio carries the denominator mean
in its own denominator, so sites with low control occupancy receive large
variances and small weights. Low occupancy is a property of the site, not noise
about the structure-response relationship, and downweighting those sites
therefore discards signal rather than noise. We suggest that reuse of
ratio-based site-level response labels treat inverse-variance weighting as a
modelling choice that requires justification, and that unweighted and weighted
analyses be reported alongside each other.

The second is condition identity, which is less deep but more easily fixed. The
archive deposits four substrate loads; the published site-level figure reports
one of them without saying which, and we identified it only by matching
replicate values, site counts and the number of highly responsive sites. A
secondary user who takes labels from that figure inherits both a condition and
the set of lysines that happened to pass the confidence threshold at that
condition, and the set differs appreciably between loads: twenty-two, thirty-
three, thirty-seven and thirty-eight sites, with only twenty passing in all
four. Stating the condition in the figure legend would cost nothing and would
remove the ambiguity entirely. More broadly, when a deposition contains several
conditions, a reuser should test across them rather than treat the published one
as representative; in this case that test was available at no experimental cost
and changed the conclusion.

With labels constructed accordingly, the biological question admits a clearer
answer than the source report could give. Local structural descriptors carried
permutation-controlled predictive information at all four substrate loads. This
is informative about the source study's own internal tension: its Results text
reported specific structural preferences, higher activity near cysteines and
reduced activity near bulky hydrophobic or aromatic residues, while its figure
legend stated that activity showed no clear preference based on primary or
secondary structure. Our result supports the former reading. The preferences are
real and detectable, they survive spatial blocking and a site-structure
permutation control, and they reproduce across conditions, which a
single-condition analysis cannot demonstrate. We emphasise that this is
quantification and confirmation of hypotheses the source authors stated, not
their discovery.

The relationship is nonetheless not uniform across sites, and we regard
establishing that as part of the contribution rather than a shortfall. The
twenty lysines quantified confidently at every load showed no
permutation-controlled signal in any condition, and in the published condition
the rank correlation was slightly negative. Random-subset resampling separates
this from a power effect: three quarters of random twenty-site subsets of the
published condition retained a rank correlation above 0.3, with a median of
0.461, so twenty sites suffice when drawn at random, and the matched set fell at
the third percentile. What distinguishes those sites we cannot say.
Measurement confidence does not predict response, the matched and unmatched
groups have indistinguishable response distributions, and no descriptor
separates them after correction for multiple testing. The one descriptor
pointing consistently in the same direction across two uncorrected tests was
neighbourhood hydrophobic fraction, which is the family the source study
associated with reduced activity; with thirty-three sites and nine descriptors
we report it as a lead and nothing more. We also note that the random subsets are
drawn from the same thirty-three sites and overlap the matched set
substantially, so the reference distribution is not independent of it, and a
third-percentile outcome arises by chance in roughly one draw in thirty. The
power explanation is excluded; the alternative, that the matched sites form a
structurally distinct class, is not established.

A second and independent instability appeared when the frozen model was compared
with the bulk responses the source study reported for casein, haemoglobin and
collagen. The experimental ordering itself is consistent, casein exceeding
haemoglobin and haemoglobin exceeding collagen at all three modification levels
deposited with the source report, so the comparison does not fail on the
measurement side. Under the prespecified specification it could not be made at
all, because no collagen lysine satisfied both the structure-confidence and
applicability-domain criteria; across the ninety-six of one hundred and
ninety-two pipeline configurations in which all three scaffolds were
represented, only seven recovered the reported ordering, against a chance
expectation of one in six. These configurations are nested variants of a single
pipeline and strongly correlated, so that fraction is descriptive rather than a
test, and selecting the agreeing minority after the fact would be specification
search. The mechanistic reading is straightforward: monomeric predicted
structures do not represent a casein micelle, the subunit interfaces of a
haemoglobin tetramer or a collagen triple helix, and a bulk immunoassay endpoint
does not identify which lysines carried the modification. We also note a hard
limit on what this comparison could ever have shown, since with three scaffolds
an exact-ordering test cannot yield a p-value below one in six.

Two boundaries define how these results should be used. All labelled sites come
from a single protein and a single assay, so spatial blocking reduces leakage
between neighbouring lysines but cannot create protein-level independence, and
nothing here licenses extrapolation to other scaffolds; the cross-scaffold
comparison above is direct evidence of that boundary rather than an argument
against it. And the four substrate loads are conditions within one experiment, so
agreement across them establishes robustness to substrate load rather than
independent biological replication, which is the sense in which we use the term
throughout.

Sample size shapes the reporting rather than the conclusions. With twenty-two to
thirty-eight sites per condition we report permutation-controlled error
differences and rank correlations rather than coefficients of determination,
which are unreliable at this scale, and we have refrained from interpreting the
ordering of conditions or from assigning causal roles to individual descriptors.
No public site-resolved response dataset exists for a second protein, so the
external comparison available to us was necessarily at protein level; that gap is
a property of the field at present rather than of this analysis, and it defines
what would most advance the question.

These limitations point to a specific next experiment rather than to a larger
model. What is needed is a prespecified panel of structurally diverse proteins,
measured at modified-lysine resolution under a shared active-versus-inactive
protocol and at more than one substrate load, with the condition reported. Such a
panel would also make it worthwhile to move beyond the single crystallographic
model used here toward conformational ensembles and native assemblies, which the
present descriptors do not represent.
Measuring more than one load is the inexpensive part and, on the evidence here,
the part that determines whether a structure-response relationship inferred
from such data can be trusted.

---

## TABLO PLANI (JPR)

- **Table 1** Data coverage: yuk, deposited rows, conf>90, used n, SHA-256 on eki
- **Table 2** Label weighting diagnostics: agirlik min/max, Kish ESS, ESS/n,
  kirpilan bolge sayisi, resp=0 ve resp=1 sayilari
- **Table 3** Condition replication: agirliksiz ve ters-varyans yan yana
- **Table 4** Matched 20-site core + rastgele altkume referansi

## SEKIL PLANI (JPR)

- **Figure 1** Reuse workflow: dort kosullu arsiv -> kosul kimligi -> etiket
  kurulumu -> kosul bazli test -> bolge kompozisyonu
- **Figure 2** Agirliklama: (A) agirlik dagilimi (log olcek, dort kosul),
  (B) Kish ESS/n, (C) agirliksiz vs ters-varyans Spearman karsilastirmasi,
  (D) permutasyon bos dagilimlari (dort kosul, agirliksiz)
- **Figure 3** Bolge kompozisyonu: (A) 200 rastgele 20-altkume Spearman
  histogrami + eslesmis cekirdegin %3 konumu, (B) cekirdekte kosul bazli
  Spearman, (C) guven-yanit sacilimi, (D) cekirdek vs disi yanit dagilimi

## YAPILACAKLAR (bu taslak sonrasi)

- [x] Yanit dagilimi 44-bolge setiyle hesaplandi ve Results'a eklendi
- [x] Tartisma bolumu yazildi (tam Limitations dahil)
- [ ] Sekiller uretilmeli
- [ ] JPR yazar kilavuzu kelime/sekil sinirlari kontrol edilmeli
- [ ] Methods (A3 taslagi) JPR cercevesine gore vurgu sirasi guncellenmeli
