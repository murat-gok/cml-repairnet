# Cok-kosullu yeniden analiz (PXD070124)

Son guncelleme: 2026-11-05. Bu klasordeki tum tablolar
`scripts/run_multicondition_reanalysis.py` tarafindan uretilir ve `--check`
bayragiyla makalede basili degerlere karsi dogrulanir.

```
python scripts/run_multicondition_reanalysis.py --check
```

Kosu kucuk bir makinede birkac bin orman egitir ve saatler surebilir. Her
tablo bitirildigi anda diske yazilir, yeniden calistirildiginda var olanlar
atlanir; yani yarida kalan bir kosu kaldigi yerden devam eder. Bir tabloyu
bastan hesaplatmak icin dosyasini silmek yeterli.

## Kaynak

`https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD070124/BSA_{0p25,0p5,1p0,2p0}_865_794_processed_filt_quant.csv`

SHA-256 on ekleri `multicondition_manifest.json` ve makalenin Tablo 1'i
icinde. Ham MS yeniden isleme YAPILMADI; depoda hazir bulunan islenmis
bolge duzeyli niceleme kullanildi.

## Etiket formulu

Yayinlanan `prepare_labels.py` ile birebir ayni: 3 aktif / 3 inaktif replike
ortalamasi, `resp = clip(1 - active_mean/control_mean, 0, 1)`, delta-metodu
varyansi `var(A)/C^2 + A^2*var(C)/C^4`, `weight = 1/(var+1e-4)` yuk icinde
medyana normalize. Bolge numaralari PRIDE "Residue #" = site_uniprot;
`site_pdb = site_uniprot - 24` (dogrulandi).

## Sadakat kontrolu

Bagimsiz yeniden uygulama, yayinlanan kosulu (500 ng) yakindan uretiyor:
MAE 0.184 (yayinlanan 0.186), Spearman 0.681 (0.699), dMAE 0.100 (0.113).
Yani sonraki farklar pipeline degisikliginden degil, verinin kendisinden
geliyor.

## Dosyalar

| Dosya | Icerik |
|---|---|
| `data_coverage.csv` | Her yukte nicelenen, guven esigini gecen ve yapisal oznitelige sahip bolge sayisi |
| `label_diagnostics.csv` | Agirlik araligi, Kish ESS, ESS/n, kirpilma sayilari (makale Tablo 2) |
| `weighting_none.csv` | Agirliksiz kosul bazli tekrar — **birincil analiz** |
| `weighting_inv_var.csv` | Ters-varyans agirliklama (yayinlanan secim) |
| `weighting_clip5.csv` | [0.2, 5] araligina kirpilmis agirliklar |
| `weighting_sensitivity.csv` | Yukaridaki uc semanin birlesimi |
| `condition_replication.csv` | Makale Tablo 3'un kaynagi (agirliksiz + ters-varyans) |
| `matched_core.csv` | Dort yukte de olculen 20 bolge, sabit kompozisyon (Tablo 4) |
| `random_subset_draws.csv` | 500 ng icinden 200 rastgele 20-bolge altkumesi |
| `response_distribution.csv` | Yuk basina yanit dagilimi |
| `confound_*.csv` | Olculebilirlik karistiricisi testleri |
| `level_labels_with_weights.csv` | Turetilmis etiket tablosu, agirliklariyla |
| `multicondition_manifest.json` | Kaynak, SHA-256'lar, calistirma parametreleri, ortak cekirdek |

`figure_data.npz` sekil betiklerinin girdisidir, bu kosu tarafindan
uretilmez.

## Ana sonuclar

**1. Agirliksiz, yapi-yanit iliskisi dort yukun dordunde de permutasyon
kontrolunu geciyor** (P = 0.010, 0.005, 0.005, 0.005). Ters-varyans
agirliklamayla 3/4; kirpilmis semayla da 3/4, ayni kosulda kaliyor. Ayrinti:
`WEIGHTING_FINDING_TR.md`.

**2. Yanit dagilimi farklari aciklamiyor.** "500 ng'de sinyal var cunku
dinamik aralik genis" hipotezi tutmuyor: sd 0.298–0.357, IQR 0.492–0.574,
yani yuklerin araliklari birbirine yakin.

**3. Dort yukte de olculen 20 bolgelik cekirdekte sinyal hicbir kosulda yok**
(P = 0.080 / 0.194 / 0.264 / 0.129). Bu bir guc kaybi degil: 500 ng icinden
cekilen 200 rastgele 20-bolge altkumesinde ortalama Spearman 0.400 (SD 0.265),
medyan 0.461, %74'u 0.3 ustunde — cekirdek bu dagilimin 3. persentilinde
(rho = -0.224). Yani sinirlayici olan orneklem buyuklugu degil, bolge
kompozisyonu.

**4. Olculebilirlik karistirici degil.** Yapisal tanimlayicilar bir bolgenin
kac kosulda guven esigini gectigini ongormuyor; cekirdek ile cekirdek disi
bolgeler yapisal olarak ayrismiyor; 500 ng icinde guven skoru yanit ile
iliskili degil.

## Model

ExtraTrees(500, max_depth=4, min_samples_leaf=2), 9 transfer-safe yapisal
oznitelik, Ward uzamsal bloklama (tam setler 4 grup, 20-bolge setleri 3
grup), on-isleme yalniz egitim katlarina fit edildi. Negatif kontrol: yapisal
oznitelik vektorlerinin bolgeler arasi permutasyonu (200 permutasyon),
etiketler ve partisyonlar sabit. 200 permutasyonda P'nin cozunurluk siniri
0.005.

## Oznitelik kapsami

Yayinlanan analizin dagittigi tablo yalnizca yayinlanan kosulun 33 bolgesini
kapsiyordu. Birlesimde yer alan 11 bolge (36, 44, 151, 160, 228, 235, 248,
285, 340, 544, 548) 4F5S'ten cikarildi ve mevcut tabloyla 33 ortak bolgede 13
alan boyunca karsilastirildi: tumu kayan nokta hassasiyetinde ayni
(max |fark| = 1.4e-14). Bu, kosullari 22/33/37/38 bolgeye cikardi ve 1000 ng
kosulunu P = 0.020'den 0.005'e, 250 ng'yi 0.030'dan 0.010'a tasidi.

Betik: `scripts/build_structure_features_union.py`
Tablo: `data/processed/bsa_structure_features_union44.csv`
