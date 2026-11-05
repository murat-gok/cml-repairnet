# Results — Proteomics cercevesine gore yeniden yazim

> **Not (2026-11-05):** Bu belge, yapisal oznitelik tablosu yalnizca
> yayinlanan kosulun 33 bolgesini kapsarken yazildi. Birlesimdeki 11 bolge
> 4F5S'ten cikarildiktan sonra kosullar 22/33/37/38 bolgeye yukseldi ve
> asagidaki sayilarin bir kismi degisti (ornegin eslesmis cekirdegin
> Spearman'i -0.131 -> -0.224, persentil %6 -> %3, ters-varyans altinda
> gecen kosul sayisi 2/4 -> 3/4). Guncel ve makaleyle uyumlu degerler icin
> `results/multicondition/` ve `docs/submission_JPR/manuscript_jpr.pdf`.
> Bu belge revizyonun izi olarak tutuluyor.

Tarih: 2026-11-03 | A1 taslaginin devami

## YENI SIRA (ve neden)

Eski Results dort alt bolumden olusuyordu ve atlas ile cross-scaffold ana
metnin yarisini aliyordu. Yeni sira tek bir anlatiya hizalandi: *yayinlanmis
kosul hangisi -> onu yeniden urettik -> digerlerinde tutmuyor -> neden
tutmadigini eledik*.

1. Kosul kimliginin saptanmasi
2. Yayinlanmis kosulun bagimsiz yeniden uretimi
3. Dort kosulda tekrar
4. Yanit dagilimlari (dinamik aralik aciklamasinin elenmesi)
5. Eslesmis bolge cekirdegi ve rastgele altkume cozumlemesi
6. Olcum karistiricisinin elenmesi

Atlas ve cross-scaffold Results'tan tamamen cikti (cross-scaffold Tartisma'da
tek paragraf, atlas supplementary).

---

# TASLAK METIN

## 3.1 The published site-level figure corresponds to the 500 ng condition

The deposited dataset contains processed site-level quantification for four
substrate loads, each with three active-enzyme and three inactive-control
replicates (Table 1). The source report does not state which load the published
site-level figure represents. Three independent lines of evidence identify it as
the 500 ng condition. First, the replicate-level percent-abundance values
underlying the published figure match that condition exactly; for example, the
lysine at position 65 carries active replicates of 9.98, 8.63 and 9.43 and
inactive replicates of 17.64, 17.40 and 17.92 in both. Second, the 500 ng
condition is the only one in which exactly thirty-three lysines pass the
confidence threshold used in the source analysis, the other loads yielding
twenty-two, thirty-seven and thirty-eight; the source report describes
"thirty-three CML modified lysine residues detected by proteomics". Third, seven
sites in that condition show more than ninety percent reduction, the number
reported in the source study.

Two counts differ slightly from the source report and we record them rather than
reconcile them. We find thirty-two of thirty-three sites with any reduction where
the source reports thirty, and twenty sites above fifty percent reduction where
the source reports twenty-one. Both discrepancies are consistent with different
handling of values at the clipping boundary and neither affects condition
identification, which rests on the exact replicate match.

Identifying the condition is a precondition for reuse rather than a result about
the enzyme. Users who take site-level labels from the published figure inherit
one substrate load together with the particular set of lysines that passed the
confidence threshold at that load, and neither inheritance is visible from the
figure itself.

## 3.2 An independent implementation reproduces the published condition

Applied to the 500 ng condition, the nine-descriptor transfer-safe
representation reached a mean absolute error of 0.184 against a mean-predictor
baseline of 0.300, with Spearman rank correlation 0.681 and a permutation-control
advantage of 0.100 (empirical P = 0.005). The corresponding values reported for
the same sites in the original analysis are 0.186, 0.699 and 0.113. The agreement
establishes that our feature definitions, spatial blocking, weighting and
permutation control reproduce the published behaviour, so the condition-level
results that follow are not attributable to implementation differences.

## 3.3 The structure-response signal is present in two of four conditions

Repeating the identical analysis in each condition separately gives a divided
result (Table 2, Figure 2). The permutation-controlled signal was recovered at
500 ng (Spearman 0.681, Delta MAE 0.100, P = 0.005) and at 2000 ng (Spearman
0.552, Delta MAE 0.086, P = 0.005). It was absent at 250 ng (Spearman 0.104,
Delta MAE 0.033, P = 0.209) and did not reach significance at 1000 ng (Spearman
0.168, Delta MAE 0.042, P = 0.060). Relative improvement over the mean predictor
followed the same pattern, ranging from 0.178 at 250 ng to 0.385 at 500 ng.

The two conditions that replicate are not adjacent in substrate load, and the
two that fail lie between and below them, so the pattern does not describe a
monotone dependence on enzyme-to-substrate ratio. With twenty-one to
thirty-three sites per condition we do not interpret the ordering of conditions.

## 3.4 Dynamic range does not explain the difference between conditions

A natural explanation for a signal appearing in some conditions only is that
response variation collapses in the others, leaving nothing to predict. The data
do not support it (Table 2, lower panel). Response standard deviations were 0.356,
0.310, 0.307 and 0.292 and interquartile ranges 0.504, 0.509, 0.430 and 0.450
across the four loads, so dispersion was comparable throughout. The conditions
differed in how responses were distributed at the upper boundary, with seven
sites above 0.90 at 500 ng against five, zero and one at the other loads, but
total spread was not reduced in the conditions where the signal failed.

## 3.5 The signal does not survive on sites measured in every condition, and this
## is not a power effect

Twenty lysines passed the confidence threshold in all four conditions and also
had extractable structural descriptors. Restricting each condition to these
twenty sites removes site composition as a difference between conditions. Under
that restriction no condition retained a permutation-controlled signal (Table 3):
empirical P values were 0.110, 0.294, 0.065 and 0.274, and in the published
condition the rank correlation became slightly negative at -0.131.

Because reducing thirty-three sites to twenty also reduces power, we tested
whether the loss is explained by sample size. Among two hundred random
twenty-site subsets of the published condition, the mean Spearman correlation was
0.414 (standard deviation 0.275), the median 0.480, and seventy-six percent of
subsets retained a correlation above 0.3; mean improvement over the baseline was
0.077. Twenty sites are therefore sufficient to recover the relationship when the
sites are drawn at random. The matched twenty-site set reached Spearman -0.131
and a baseline improvement of 0.009, placing it at the sixth percentile of the
random-subset distribution (Figure 3A).

A sixth-percentile outcome occurs by chance in roughly one draw in seventeen, so
this analysis establishes that power is not the explanation without establishing
that the matched sites are a distinct class. The signal in the published
condition is carried disproportionately by the thirteen lysines that are not
confidently quantified at every substrate load.

## 3.6 Measurement confidence does not confound the relationship

Because the thirteen sites that carry the signal are those quantified less
consistently across loads, we tested whether structural descriptors predict
measurability rather than response. They do not, after correction for testing
nine descriptors (Bonferroni threshold 0.0056). The strongest association with
the number of conditions in which a site was confidently quantified was
neighbourhood hydrophobic fraction (Spearman 0.442, P = 0.0101), which does not
survive correction, and no descriptor was associated with mean confidence (all
P > 0.11). Within the published condition, confidence score and response were
unrelated (Spearman -0.087, P = 0.63). The matched and unmatched site groups had
indistinguishable response distributions, with means of 0.602 and 0.603
(Mann-Whitney P = 0.84).

Comparing the two site groups descriptor by descriptor gave no difference
surviving correction either, the smallest P value again belonging to
neighbourhood hydrophobic fraction, which averaged 0.404 in matched and 0.251 in
unmatched sites (P = 0.0164), followed by distance to the protein centroid
(29.17 and 24.36, P = 0.068) and residue solvent-accessible surface area (88.21
and 67.86, P = 0.126). We report the convergence of hydrophobic neighbourhood
across two uncorrected tests as a lead, noting that it is the descriptor family
the source study associated with reduced activity, and not as a finding.

---

## TABLOLAR

### Table 1. Data coverage across the four substrate loads

| Substrate load (ng) | File key | Rows deposited | Sites above confidence threshold | Sites with structural descriptors | Used |
|---|---|---|---|---|---|
| 250 | 0p25 | 35 | 22 | 21 | 21 |
| 500 | 0p5 | 40 | 33 | 33 | 33 |
| 1000 | 1p0 | 44 | 37 | 29 | 29 |
| 2000 | 2p0 | 45 | 38 | 30 | 30 |

Footnote: union across loads 44 sites; 20 sites pass in all four loads;
structural descriptors available for 33 sites; 11 sites of the union
(36, 44, 151, 160, 228, 235, 248, 285, 340, 544, 548) lacked extractable
descriptors under the quality rules and were excluded. Each condition has three
active and three inactive replicates. SHA-256 prefixes of the four source files
are given in the data availability statement.

### Table 2. Condition-level replication and response distribution

| Load (ng) | n | MAE | Baseline MAE | Relative improvement | Spearman | Permuted MAE | Delta MAE | P |
|---|---|---|---|---|---|---|---|---|
| 250 | 21 | 0.327 | 0.398 | 0.178 | 0.104 | 0.360 | 0.033 | 0.209 |
| 500 | 33 | 0.184 | 0.300 | 0.385 | 0.681 | 0.284 | 0.100 | 0.005 |
| 1000 | 29 | 0.240 | 0.307 | 0.220 | 0.168 | 0.282 | 0.042 | 0.060 |
| 2000 | 30 | 0.200 | 0.280 | 0.287 | 0.552 | 0.286 | 0.086 | 0.005 |

| Load (ng) | n | Mean response | SD | Min | Max | IQR | Sites <= 0.10 | Sites >= 0.90 |
|---|---|---|---|---|---|---|---|---|
| 250 | 21 | 0.530 | 0.356 | 0.00 | 0.96 | 0.504 | 5 | 5 |
| 500 | 33 | 0.602 | 0.310 | 0.00 | 1.00 | 0.509 | 3 | 7 |
| 1000 | 29 | 0.504 | 0.307 | 0.00 | 0.89 | 0.430 | 5 | 0 |
| 2000 | 30 | 0.489 | 0.292 | 0.00 | 0.97 | 0.450 | 3 | 1 |

### Table 3. Matched twenty-site core and random-subset reference

| Load (ng) | n | MAE | Baseline MAE | Spearman | Delta MAE | P |
|---|---|---|---|---|---|---|
| 250 | 20 | 0.313 | 0.511 | -0.035 | 0.051 | 0.110 |
| 500 | 20 | 0.255 | 0.264 | -0.131 | 0.019 | 0.294 |
| 1000 | 20 | 0.214 | 0.246 | 0.062 | 0.031 | 0.065 |
| 2000 | 20 | 0.259 | 0.317 | 0.124 | 0.020 | 0.274 |

Random-subset reference, 200 random twenty-site subsets of the 500 ng condition:
Spearman mean 0.414, SD 0.275, 5th percentile -0.159, median 0.480, 95th
percentile 0.724; 76% of subsets above 0.3; mean baseline improvement 0.077.
Matched core: Spearman -0.131, baseline improvement 0.009, 6th percentile.

---

## SEKIL PLANI

Mevcut Sekil 2 ve 3 bu cerceveye uymuyor (tek kosullu). Yeni plan:

**Figure 1 — Reuse workflow.** Dort kosullu arsivden kosul kimligi saptamasina,
oradan kosul-bazli tekrar ve bolge-kompozisyonu cozumlemesine akis. Eski
Sekil 1 yeniden cizilecek; atlas kolu cikacak.

**Figure 2 — Condition-level replication.** (A) Her kosul icin gozlenen-tahmin
sacilimi, 2x2 panel. (B) Her kosul icin permutasyon bos dagilimi (histogram) ve
gercek MAE'nin konumu — dort kucuk panel; 0p5 ve 2p0'da gercek deger dagilimin
solunda, 0p25 ve 1p0'da icinde. (C) Kosul basina yanit dagilimi (violin veya
kutu). (D) Kosul basina Spearman, bootstrap aralikli.

**Figure 3 — Site composition, not power.** (A) 200 rastgele 20-bolge
altkumesinin Spearman histogrami, eslesmis cekirdegin konumu dikey cizgi ile
(%6 percentile) — makalenin en onemli tek paneli. (B) Eslesmis cekirdekte kosul
basina Spearman. (C) 0p5 icinde guven skoru ve yanit sacilimi (iliskisizlik).
(D) Eslesmis vs eslesmeyen bolgelerde yanit dagilimi (ayirt edilemezlik) ve
hidrofobik komsuluk orani (tek yonlu iz).

**Supplementary.** Eski Sekil 4 (transfer limitleri) ve atlas sekilleri buraya
iner. Eski Sekil 2-3 panellerinin 0p5'e ait olanlari Figure 2'ye devrolur.

---

## TUTARLILIK KONTROL LISTESI (metin yazilirken)

- [ ] "framework" kelimesi hicbir yerde yok
- [ ] R-kare hicbir yerde yok
- [ ] Kosullar "substrate load" olarak adlandirildi, "enzyme:substrate ratio" degil
- [ ] Her P degeri permutasyon tabanli ve kac permutasyon oldugu belirtilmis
- [ ] Ters-varyans agirliklama Methods'ta aciklanmis (B1 maddesi)
- [ ] Kish etkin orneklem buyuklugu raporlanmis (B2 maddesi)
- [ ] Atlas ana metinde gecmiyor
- [ ] 11 eksik bolge ve nedeni Table 1 dipnotunda
- [ ] Hidrofobik iz "lead, not finding" olarak nitelenmis
- [ ] Iki kucuk sayim tutarsizligi (30 vs 32, 21 vs 20) aciklanmis
