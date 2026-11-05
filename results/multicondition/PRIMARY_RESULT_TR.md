# BIRINCIL SONUC — 44-bolge oznitelik seti, agirliksiz

Son guncelleme: 2026-11-05. Bu, makalenin ana bulgusunun nihai hali.

| substrat yuku (ng) | n | MAE | baseline | rel. iyilesme | Spearman | dMAE | P |
|---|---|---|---|---|---|---|---|
| 250  | 22 | 0.274 | 0.337 | 0.189 | 0.414 | 0.090 | **0.010** |
| 500  | 33 | 0.187 | 0.275 | 0.319 | 0.616 | 0.102 | **0.005** |
| 1000 | 37 | 0.221 | 0.268 | 0.176 | 0.399 | 0.064 | **0.005** |
| 2000 | 38 | 0.188 | 0.270 | 0.303 | 0.581 | 0.091 | **0.005** |

**Dort kosulun dordunde de permutasyon kontrolu geciliyor** (P <= 0.010; ucu
200-permutasyon tabani olan 0.005'te).

## Eksik bolgelerin eklenmesinin etkisi

| yuk | eski n | eski P | yeni n | yeni P |
|---|---|---|---|---|
| 250  | 21 | 0.0299 | 22 | **0.0100** |
| 500  | 33 | 0.0050 | 33 | 0.0050 |
| 1000 | 29 | 0.0199 | 37 | **0.0050** |
| 2000 | 30 | 0.0050 | 38 | 0.0050 |

11 eksik bolgenin 4F5S'ten cikarilmasi 1000 ng kosulunu P=0.020'den 0.005'e,
250 ng'yi 0.030'dan 0.010'a tasidi. Bedeli sifir: PDB zaten repoda idi.

## Oznitelik cikariminin dogrulanmasi

Cikarim mevcut `bsa_structure_features.csv` ile 33 ortak bolgede 13 oznitelik
icin karsilastirildi: tum alanlar kayan nokta hassasiyetinde ayni
(max |fark| = 1.4e-14). Betik: `scripts/build_structure_features_union.py`
Tablo: `data/processed/bsa_structure_features_union44.csv`

## Duyarlilik ve kontroller (tamamlandi)

- **Ters-varyans, 44-bolge seti:** 3/4 kosul geciyor. 250 ng kaliyor
  (Spearman 0.141, P = 0.109); 500/1000/2000 geciyor (P = 0.005 / 0.010 /
  0.005).
- **[0.2, 5] kirpilmis:** yine 3/4, ayni kosulda kaliyor — kirpma agirlikli
  sonucu yeniden uretiyor, agirliksizi degil.
- **Eslesmis 20-bolge cekirdegi, agirliksiz:** hicbir kosulda sinyal yok
  (P = 0.080 / 0.194 / 0.264 / 0.129).
- **Rastgele 20-altkume referansi (500 ng):** ortalama Spearman 0.400
  (SD 0.265), medyan 0.461, %74'u 0.3 ustunde. Cekirdek 3. persentilde
  (rho = -0.224), yani cekirdek sonucu guc kaybi degil, kompozisyon.
- **Olculebilirlik karistiricisi:** yapisal tanimlayicilar bir bolgenin kac
  kosulda esigi gectigini ongormuyor; cekirdek ile disi yapisal olarak
  ayrismiyor.

## Yeniden uretim

```
python scripts/run_multicondition_reanalysis.py --check
```
