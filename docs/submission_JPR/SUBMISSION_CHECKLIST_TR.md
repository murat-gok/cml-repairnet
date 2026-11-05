# JPR gonderim kontrol listesi — FINAL
Guncelleme: 2026-11-05

## TAMAMLANDI

- [x] `manuscript_jpr.tex` / `.pdf` — achemso (`journal=jprobs`), 26 sayfa
- [x] `supporting_information.tex` / `.pdf` — 3 sayfa (S1 protein duzeyi
      karsilastirma + Figure S1, S2 yanit dagilimi, S3 tanimlayici dogrulamasi)
- [x] `references_jpr.bib` — 22 girdi, **tum yazar listeleri tam**
- [x] 3 ana sekil + TOC graphic + Figure S1: vektor PDF + 600 dpi LZW-TIFF
- [x] 4 tablo (Tablo 3 her iki agirliklama semasi icin MAE/rho/P veriyor)
- [x] TOC graphic 3.25 x 1.75 inc (achemso kutusuyla birebir)
- [x] Keywords (6), Author Information, Author Contributions, Notes (COI),
      Data and Code Availability

## SON KONTROLDE BULUNAN VE DUZELTILEN

1. **Methods 2.3 eski 33-bolge sayilarini tasiyordu** — "0.24 ve 0.31" ve
   "two of four conditions" yaziyordu, oysa Ozet dogru olarak "0.23-0.32" ve
   "three of four" diyordu. Abstract ile Methods celisiyordu. Duzeltildi.
2. **Sekiller metinde hic atif almiyordu** — Figure 1/2/3 icin bes atif eklendi.
3. **Tablo 3'te ham MAE yoktu** — yalnizca dMAE vardi; MAE ve baseline eklendi.
4. **Alti kaynakta yazar listesi kisaltilmisti** — Crossref API ve JMLR resmi
   sayfasindan dogrulanip tamamlandi (PRIDE 11 yazar, scikit-learn 16,
   SciPy 34 + konsorsiyum, Roberts 14, Chaudhuri 8, Biopython 11).

## OTOMATIK DOGRULAMA (son derleme)

- Cozulmemis referans (`??`) : 0
- Cozulmemis atif (`(?)`)    : 0
- Kaynakcada "others"        : 0
- Eski 33-bolge degerleri    : 0 (tarandi: -0.131, sixth percentile, 0.0299,
                               0.0199, 0.446, 0.436, twenty-one to thirty-three)
- Yasakli kalintilar         : 0 ("framework", R-kare, atlas, Pareto, STRING,
                               g:Profiler, "enzyme:substrate ratio")
- Tum sekil ve tablolar metinde atif aliyor

## GONDERMEDEN ONCE YAZARIN YAPMASI GEREKEN

ACS'nin resmi sekil/LaTeX kural sayfalarina erisilemedi (404 / Cloudflare 403).
Sekiller en kati yaygin gereklilige gore uretildi: vektor PDF + 600 dpi TIFF,
tek sutun 8.3 cm / cift sutun 17.8 cm, Liberation Sans (Arial ile metrik
uyumlu), en kucuk etiket 6 pt, RGB, gomulu font. Submission sistemi farkli
sayilar isterse `scripts/make_jpr_figures.py` ve `make_jpr_toc.py` basindaki
W1/W2/DPI_TIFF/FONT degistirilip yeniden calistirilir.

## YENIDEN URETIM

```
python scripts/build_jpr_tables.py
python scripts/build_jpr_tex.py
sh     scripts/fix_jpr_tex.sh manuscript_jpr.tex
python scripts/add_jpr_citations.py
latexmk -pdf manuscript_jpr.tex
latexmk -pdf supporting_information.tex
```
Kaynak metin `docs/submission_JPR/MANUSCRIPT_JPR_TR.md` ve
`docs/revision_A1/A3_METHODS_PROTEOMICS_TR.md`; `.tex` elle duzenlenmemeli.
