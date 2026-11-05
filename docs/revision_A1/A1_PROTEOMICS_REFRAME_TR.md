# A1 — Proteomics (Wiley) icin yeniden cerceveleme taslagi

> **Not (2026-11-05):** Bu belge, yapisal oznitelik tablosu yalnizca
> yayinlanan kosulun 33 bolgesini kapsarken yazildi. Birlesimdeki 11 bolge
> 4F5S'ten cikarildiktan sonra kosullar 22/33/37/38 bolgeye yukseldi ve
> asagidaki sayilarin bir kismi degisti (ornegin eslesmis cekirdegin
> Spearman'i -0.131 -> -0.224, persentil %6 -> %3, ters-varyans altinda
> gecen kosul sayisi 2/4 -> 3/4). Guncel ve makaleyle uyumlu degerler icin
> `results/multicondition/` ve `docs/submission_JPR/manuscript_jpr.pdf`.
> Bu belge revizyonun izi olarak tutuluyor.

Tarih: 2026-11-03

## ONCE IKI DUZELTME

**1. Kosullar oran degil, substrat MIKTARI.** Trabosh Methods: "Either 250, 500,
1000, or 2000 ng CML-modified BSA in 25 uL reaction volume was incubated with
enzyme samples at 5 uM concentration for 20 h at 24 C." Yani enzim sabit,
substrat artiyor. Dosya adlari 0p25/0p5/1p0/2p0 = 0.25/0.5/1.0/2.0 ug substrat.
Artan substrat = DUSEN enzim:substrat orani. Onceki yazismalarda bunlara
"enzim:substrat orani" dedim — yanlisti, metinde "substrate load" kullanilacak.

**2. Yayinlanmis sekil 500 ng kosuluna ait.** Makale hangi kosulu gosterdigini
BELIRTMIYOR ("Neither the Figure 3B nor the Figure 3C legend specifies which
substrate amount"). Uc bagimsiz kanitla 0p5 oldugu saptandi:
- `bsa_labels.csv` replike degerleri 0p5 dosyasiyla birebir ayni (orn. bolge 65:
  aktif 9.98/8.63/9.43, inaktif 17.64/17.40/17.92)
- 0p5 conf>90 bolge sayisi = 33; makale "thirty-three CML modified lysine
  residues detected" diyor; diger kosullar 22/37/38
- 0p5'te >%90 azalma gosteren bolge sayisi = 7; makale 7 bildiriyor
Bu, makalede olmayan ve yeniden kullanim icin gerekli bir bilgi.

## TRABOSH'UN KENDI NITELIKSEL IDDIALARI (birebir alinti, A1'in dayanagi)

Bunlar Introduction'da ACIKCA atif alacak — makalenin merkezi sorusu bu
gozlemlerin nicelleştirilmesi oldugu icin:

- "CMLase activity on a particular site was not strictly a function of surface
  accessibility"
- K51, K232, K573 icin: "all reside in surface-exposed helices but showed varied
  CML reductions of 90%, 38%, and 0%"
- "Activity weakly correlated with local protein flexibility (B-factor)" /
  "flexibility may be necessary but is not sufficient for catalysis"
  (Supplementary Fig. 13)
- "Some local sequence preference was observed, with CMLase showing higher
  activity on sites with neighboring cysteines" (Supplementary Fig. 14)
- "Reduced enzymatic activity was associated with bulky hydrophobic or aromatic
  residues" (Supplementary Fig. 14)
- Sekil 3C aciklamasi: "CMLase has varying activity across all sites with no
  clear preference based on primary or secondary structure."

**Kullanilacak ic gerilim:** Results metni belirli tercihler bildiriyor
(Cys+ , hacimli hidrofobik/aromatik -), Sekil 3C aciklamasi ise "no clear
preference" diyor. Bizim nicel sonucumuz bu gerilimi cozuyor: tercih VAR ama
substrat yukune gore istikrarsiz. Bu, makalenin en savunulabilir katkisi.

---

# TASLAK METIN

## Baslik

**Site-level carboxymethyllysine repair response is not stable across substrate
load: a four-condition reanalysis of a public proteomics dataset**

(Alternatif, daha yontem-vurgulu: "Reuse of single-condition site-resolved PTM
response data overstates structure-response relationships: a worked example on
enzymatic carboxymethyllysine repair")

## Ozet (~210 kelime)

Nepsilon-carboxymethyllysine (CML) is a stable advanced glycation end product
that was long considered irreversible until the recent description of CMLase, an
engineered oxidase that restores native lysine. The source study reported
site-resolved CML reduction across thirty-three bovine serum albumin lysines and
noted qualitative structural trends, including higher activity near cysteines and
lower activity near bulky hydrophobic residues, while also stating that activity
showed no clear structural preference. We reused the deposited dataset
(ProteomeXchange PXD070124) to test whether local structural descriptors carry
reproducible information about this site-specific response. The archive contains
four substrate-load conditions; the published site-level figure corresponds to
the 500 ng condition, which the source report does not specify and which we
identify here from replicate-level agreement, site count and the number of highly
responsive sites. Using a nine-descriptor transfer-safe representation with
spatially blocked cross-validation and a site-structure permutation control, we
reproduced the published condition closely (mean absolute error 0.184; Spearman
0.681) and recovered a permutation-controlled signal in two of four conditions.
The signal was absent in the remaining two and did not survive on the twenty
sites measured at high confidence in all four conditions. Random-subset analysis
excluded power loss as the explanation, and measurement confidence did not
predict response. A protein-level comparison against bulk responses for casein, haemoglobin and
collagen also failed to recover the reported ordering, and we show that with
three scaffolds such a comparison cannot reach significance in principle.
Site-level structure-response relationships inferred from a
single published condition are therefore unstable, and reuse of such data
requires condition-level replication before transferable models are claimed.

## Giris

Advanced glycation end products accumulate when proteins react non-enzymatically
with reducing sugars and reactive carbonyl species, and among them
Nepsilon-carboxymethyllysine (CML) is chemically stable, abundant in long-lived
proteins and capable of altering lysine charge and receptor-mediated signalling.
Until recently CML was treated as a terminal modification: endogenous defences
such as the glyoxalase system and fructosamine-3-kinase act on precursors or
early glycation adducts rather than on stable CML, and the reported deglycase
activity of DJ-1/PARK7 remains contested. The description of CMLase, an engineered
FAD-dependent oxidase derived from a glycine oxidase scaffold that oxidises the
Nepsilon-carboxymethyl group and restores native lysine, therefore created the
first opportunity to ask which modified lysines an enzyme can actually reach.

That source study answered the question in bulk and began to answer it at site
level. Of thirty-three CML-modified lysines detected on bovine serum albumin,
thirty showed reduced CML, but the responses spanned almost the full scale, and
the authors reported specific qualitative trends: activity was "not strictly a
function of surface accessibility", three lysines in surface-exposed helices
"showed varied CML reductions of 90%, 38%, and 0%", activity "weakly correlated
with local protein flexibility (B-factor)" such that "flexibility may be
necessary but is not sufficient for catalysis" (their Supplementary Fig. 13), and
"some local sequence preference was observed, with CMLase showing higher activity
on sites with neighboring cysteines" while "reduced enzymatic activity was
associated with bulky hydrophobic or aromatic residues" (their Supplementary
Fig. 14). The legend of their Fig. 3C nonetheless states that "CMLase has varying
activity across all sites with no clear preference based on primary or secondary
structure." Our starting point is this unresolved tension within the source
report, not a gap in it: the qualitative structural observations have been made,
and what is missing is a quantitative, leakage-controlled test of whether they
hold.

Testing them is a data-reuse problem as much as a modelling problem. The
deposited dataset contains processed site-level quantification for four substrate
loads (250, 500, 1000 and 2000 ng of CML-modified BSA incubated with enzyme at
fixed concentration), so the enzyme-to-substrate ratio falls as load rises. The
published site-level figure reports one of these conditions without identifying
which, and secondary users who take the figure source data inherit that single
condition together with its particular set of confidently quantified sites. We
therefore treated condition identity as something to establish rather than
assume, and condition-level replication as the primary test rather than an
afterthought.

This matters beyond one enzyme. Site-resolved post-translational-modification
response measurements are increasingly used as labels for structure-based
prediction, and the available label sets are small, single-protein and
single-condition. Established predictors of lysine glycation occurrence, from
early neural-network models through recent transformer-based classifiers, predict
which lysines become modified; none predicts how a modified lysine responds to an
enzyme, so there is no precedent for how stable such response labels are. If a
structure-response relationship fitted to one published condition does not hold
in sibling conditions from the same experiment, that is a property of the data
class, not of this dataset alone.

Here we reuse PXD070124 across all four substrate loads. We identify which
condition the published figure reports, reproduce the published condition with an
independent implementation, and then test a compact nine-descriptor structural
representation under spatially blocked cross-validation with a site-structure
permutation control in each condition separately. We further test the twenty
sites quantified at high confidence in all four conditions, use random-subset
resampling to separate statistical power from site composition, and test whether
measurement confidence confounds the apparent relationship. The result is a
bounded negative: local structure carries permutation-controlled information in
some conditions and not others, and the data cannot determine why.

## Tartisma

Our reimplementation reproduced the published condition closely. On the 500 ng
condition the frozen model reached a mean absolute error of 0.184 against a
mean-predictor baseline of 0.300, with Spearman 0.681 and a permutation-control
advantage of 0.100, compared with 0.186, 0.699 and 0.113 in the original analysis
of the same sites. This agreement matters methodologically rather than
scientifically: it establishes that the negative results that follow are not
artefacts of a different pipeline, feature definition or validation scheme.

Across conditions the relationship was unstable. The permutation-controlled
signal was present at 500 ng (P = 0.005) and at 2000 ng (Spearman 0.552,
P = 0.005), absent at 250 ng (Spearman 0.104, P = 0.209) and did not reach
significance at 1000 ng (P = 0.060). Two of four conditions therefore replicate
the published relationship and two do not. We considered and rejected the most
obvious explanation: response distributions were comparable across conditions
(standard deviations 0.292 to 0.356, interquartile ranges 0.430 to 0.509), so the
conditions in which the signal failed were not conditions in which dynamic range
had collapsed.

The matched-site analysis is harder to explain away. On the twenty lysines
quantified at high confidence in all four conditions, no condition retained a
permutation-controlled signal, and in the published condition the rank
correlation became slightly negative. Random-subset resampling excluded loss of
statistical power: among 200 random twenty-site subsets of the published
condition, 76% retained Spearman above 0.3 and the median was 0.480, so twenty
sites suffice to detect the relationship when the sites are drawn at random. The
matched subset sits at the sixth percentile of that distribution. We note,
however, that a sixth-percentile outcome arises by chance roughly once in
seventeen draws, so the matched subset may itself be an unusual sample rather
than a structurally distinct class of site.

We tested, and did not find, a measurement-based explanation. Structural
descriptors did not predict in how many conditions a site was confidently
quantified after correction for multiple testing, measurement confidence did not
predict response within the published condition (Spearman -0.087, P = 0.63), and
the matched and non-matched site groups had indistinguishable response
distributions (means 0.602 and 0.603, P = 0.84). The apparent relationship is
therefore not an artefact of structure determining detectability.

One descriptor pointed consistently in the same direction without surviving
correction. Neighbourhood hydrophobic fraction was higher in matched than
non-matched sites (0.404 versus 0.251) and correlated with the number of
conditions in which a site was confidently quantified. This is the descriptor
family that the source study associated with reduced activity, so the convergence
is suggestive, but with thirty-three sites and nine descriptors neither test
survives Bonferroni correction and we report it as a lead rather than a finding.

A second and independent instability appeared when the frozen model was compared
with the bulk responses the source study reported for casein, haemoglobin and
collagen. The experimental ordering itself is consistent: casein exceeded
haemoglobin, which exceeded collagen, at all three modification levels deposited
with the source report, so the comparison does not fail on the measurement side.
It fails on the prediction side. Under the prespecified primary specification the
comparison could not be made at all, because no collagen lysine satisfied both
the structure-confidence and applicability-domain criteria. Across the 96 of 192
pipeline configurations in which all three scaffolds were represented, only
seven, or 7.3%, recovered the reported ordering, against a chance expectation of
one in six. Because these configurations are nested variants of a single pipeline
they are strongly correlated, so this fraction is descriptive rather than a
significance test, and selecting the agreeing minority after the fact would be
specification search. The mechanistic reading is straightforward: monomeric
predicted structures do not represent a casein micelle, the subunit interfaces of
a haemoglobin tetramer or a collagen triple helix, and a bulk immunoassay
endpoint does not identify which lysines carried the modification. We also note a
hard limit on what this comparison could ever have shown. With three scaffolds an
exact-ordering test cannot yield a p-value below one in six, so protein-level
agreement of this kind cannot support a transfer claim however it turns out.

Two implications follow for reuse of this data class. First, a published
site-level figure is a single realisation of an experiment that may contain
several, and the condition it reports should be stated by the original authors or
established by the reuser; here it had to be inferred. Second, and more
generally, a structure-response relationship fitted to one condition should not be
described as reproducible until it has been tested in the sibling conditions that
the archive already contains. In this dataset that test was available at no
experimental cost and changed the conclusion.

The study has clear limits. All labelled sites come from one protein and one
assay, so spatial blocking reduces within-protein leakage but cannot create
protein-level independence. Twenty to thirty-three sites per condition leave the
condition-dependence and sampling-variability explanations unresolved, and we do
not claim to distinguish them. Eleven sites in the union across conditions lacked
extractable structural descriptors under our quality rules and were excluded. The
endpoint is a relative reduction in CML signal and may combine catalytic
susceptibility with modification stoichiometry and measurement effects. The most
informative next experiment is therefore not a larger model but a prespecified
panel of structurally diverse proteins measured at modified-lysine resolution
under a shared active-versus-inactive protocol and at more than one substrate
load, so that condition dependence can be estimated rather than assumed.

## Methods'a eklenecek paragraf (kisa)

Protein-level comparison. The frozen model was compared with the bulk
CML-ELISA responses reported for casein, haemoglobin and collagen in the source
study, using canonical bovine components with stoichiometric weighting because
the commercial species and isoforms were not specified. The observed side was
read at full resolution: the source report deposits three modification levels per
scaffold, and the ordering was evaluated within each level rather than after
collapsing levels to a scaffold median, since a three-point rank correlation
cannot attain a p-value below one in six. Predicted scaffold scores are
level-invariant, so the independent units for the model-observation comparison
are the three scaffolds, and the attainable p-value floor for an exact-ordering
match is therefore one in six regardless of how many levels are measured. The
comparison was prespecified as an applicability stress test, not as site-resolved
external validation, because the endpoint is a bulk immunoassay that does not
identify which lysines carried CML.

---

## NE CIKARILIYOR

- "framework", "we develop a framework" dili tamamen gitti
- R-kare hicbir yerde gecmiyor (ozet ve ana metin MAE + Spearman + permutasyon P)
- 1.268 bolgelik insan/fare atlasi ANA METINDEN cikti; Proteomics kapsamina
  girmiyor ve savunulamiyor (739 delivery-unresolved, insan enrichment'i yok,
  yuksek baglamli 10 bolgenin 7'si fare albumini = egitim iskeleti)
- Pareto/NSGA-II, STRING, g:Profiler katmanlari cikti
- "prioritizing age-related protein damage" iddiasi cikti
- Cross-scaffold (casein/Hb/kollajen) bolumu TEK PARAGRAFA indirildi
  (Tartisma'da, "Two implications follow" paragrafindan hemen once) + Methods'a
  kisa bir yontem paragrafi. Ayri Results alt bolumu, Tablo ve Sekil 4 ana
  metinden cikiyor; Sekil 4 supplementary'ye iner.

## NE EKLENIYOR

- Kosul kimligi saptamasi (yayinlanmis seklin hangi kosula ait oldugu)
- Dort kosullu tekrar testi
- Eslesmis 20-bolge testi + rastgele altkume cozumlemesi
- Olcum karistiricisi testi (negatif)
- Veri yeniden kullanimi cercevesi (Proteomics kapsamina dogrudan hitap eder)
