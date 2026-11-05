# Agirliklama bulgusu

Son guncelleme: 2026-11-05 — tum sayilar 44-bolge oznitelik setiyle, 200
permutasyonla yeniden uretildi ve `scripts/run_multicondition_reanalysis.py
--check` ile makaleye karsi dogrulandi. (Bu notun 2026-11-04 tarihli ilk hali
33-bolgelik ara sonuclari tasiyordu; asagidaki tablolar onlarin yerine gecer.)

## Sonuc

| yuk (ng) | n | inv_var rho | inv_var P | **none rho** | **none P** | clip5 rho | clip5 P |
|---|---|---|---|---|---|---|---|
| 250  | 22 | 0.141 | 0.109 | **0.414** | **0.010** | 0.299 | 0.060 |
| 500  | 33 | 0.681 | 0.005 | **0.616** | **0.005** | 0.663 | 0.005 |
| 1000 | 37 | 0.298 | 0.010 | **0.399** | **0.005** | 0.290 | 0.020 |
| 2000 | 38 | 0.574 | 0.005 | **0.581** | **0.005** | 0.563 | 0.005 |

- **Agirliksiz: 4/4 kosulda permutasyon kontrolu geciliyor.**
- Ters-varyans: 3/4 — 250 ng'de kaliyor
- [0.2, 5] araligina kirpilmis: 3/4 — ayni kosulda kaliyor, yani kirpma
  agirlikli sonucu yeniden uretiyor, agirliksizi degil

200 permutasyonda P'nin cozunurluk siniri 0.005.

## Yorum

Kosullar arasi istikrarsizlik buyuk olcude agirliklama semasinin urunu,
biyolojik degil.

Mekanizma: delta-metodu varyansi `var(A)/C^2 + A^2*var(C)/C^4` paydasinda
kontrol ortalamasi (C) tasiyor. Dusuk kontrol bollugu olan bolgeler cok buyuk
varyans, dolayisiyla cok kucuk agirlik aliyor. Ama dusuk bolluk, yapi-yanit
iliskisi hakkinda gurultu DEGIL — bolgenin bir ozelligi. Yani ters-varyans
agirliklama dusuk-doluluklu bolgeleri sistematik olarak bastiriyor ve doluluga
bagli bir yanlilik getiriyor.

Kirpmanin kurtarmamasi onemli: sorun yalnizca asiri agirliklar degil,
agirliklamanin YONU. Orta siddette ters-varyans agirliklama da zarar veriyor.
Kirpilmis sema, basarisiz olan kosulu ters-varyansla ayni birakiyor.

## Etkin orneklem (Kish)

| yuk | n | agirlik araligi | Kish ESS | ESS/n |
|---|---|---|---|---|
| 250  | 22 | 0.0002–28.4 | 5.0  | 0.23 |
| 500  | 33 | 0.0048–18.7 | 9.7  | 0.30 |
| 1000 | 37 | 0.0016–18.6 | 11.8 | 0.32 |
| 2000 | 38 | 0.0171–38.0 | 11.4 | 0.30 |

Ters-varyans agirliklama etkin orneklemi nominalin %23–32'sine indiriyor;
yayinlanan "33 bolge" kosulunda etkin orneklem 9.7. Agirlik araligi en genis
yukte bes mertebeye yayiliyor.

## Karar

**Birincil analiz AGIRLIKSIZ.** Ters-varyans ve kirpilmis semalar duyarlilik
analizi olarak raporlanir, ve dusuk-doluluk yanliligi mekanizmasi aciklanir.
Bu, makale icin hem daha saglam bir sonuc hem de yayimlanabilir bir
metodolojik gozlem: oran-temelli PTM yanit etiketlerinde ters-varyans
agirliklama, dusuk-doluluklu bolgeleri bastirarak gercek bir yapi-yanit
sinyalini yok edebilir.

## Agirliksiz olarak tekrarlananlar

Ilk sonuclarin tamami ters-varyans agirliklariyla uretilmisti; hepsi
agirliksiz olarak yeniden calistirildi:

- [x] Kosul bazli tekrar (yukaridaki tablo)
- [x] Eslesmis 20-bolge cekirdegi, dort kosul — hicbirinde sinyal yok
      (P = 0.080 / 0.194 / 0.264 / 0.129)
- [x] Rastgele 20-altkume referansi (500 ng) — ortalama rho 0.400 (SD 0.265),
      medyan 0.461, %74'u 0.3 ustunde; cekirdek 3. persentilde (rho = -0.224).
      Yani cekirdek sonucu guc kaybi degil, bolge kompozisyonu.
- [x] Karistirici testleri — yapisal oznitelikler olculebilirligi ongormuyor
      (en kucuk p, dokuz tanimlayici icinde Bonferroni esiginin uzerinde)

## Yeniden uretim

```
python scripts/run_multicondition_reanalysis.py --check
```

Tablolar `results/multicondition/` altina yazilir; `--check` bayragi uretilen
degerleri makalede basili olanlarla karsilastirir ve biri kaymissa sifirdan
farkli kodla cikar.
