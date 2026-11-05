# NIHAI SONUC TABLOSU — 44-bolge oznitelik seti (2026-11-04)

Tum analizler tamamlandi. Bu, makalede raporlanacak nihai sayilar.

## 1. Birincil analiz: agirliksiz, kosul bazli tekrar

| yuk (ng) | n | MAE | baseline | rel. iyilesme | Spearman | dMAE | P |
|---|---|---|---|---|---|---|---|
| 250  | 22 | 0.274 | 0.337 | 0.189 | 0.414 | 0.090 | **0.010** |
| 500  | 33 | 0.187 | 0.275 | 0.319 | 0.616 | 0.102 | **0.005** |
| 1000 | 37 | 0.221 | 0.268 | 0.176 | 0.399 | 0.064 | **0.005** |
| 2000 | 38 | 0.188 | 0.270 | 0.303 | 0.581 | 0.091 | **0.005** |

**4/4 kosulda permutasyon kontrolu geciliyor** (P <= 0.010; ucu 200-permutasyon
tabani olan 0.005'te).

## 2. Duyarlilik: ters-varyans agirliklama

| yuk (ng) | n | Spearman | dMAE | P | geciyor |
|---|---|---|---|---|---|
| 250  | 22 | 0.141 | 0.048 | 0.1095 | HAYIR |
| 500  | 33 | 0.681 | 0.100 | 0.0050 | evet |
| 1000 | 37 | 0.298 | 0.047 | 0.0100 | evet |
| 2000 | 38 | 0.574 | 0.079 | 0.0050 | evet |

Ters-varyans ile 3/4. Eksik bolgelerin eklenmesi agirlikli analizi de
iyilestirdi (1000 ng: P 0.060 -> 0.010), ama agirliksiz hala UNIFORM olarak
daha iyi: Spearman dort kosulun ucunde yuksek, ve tek basarisizlik en kucuk
kosulda.

Agirliklama bulgusu gecerliligini koruyor ama YUMUSADI: onceki 33-bolge
setinde 2/4 idi, simdi 3/4. Metinde bu duzeltilmis haliyle yazilacak —
"agirliklama sinyali iki kosulda yok ediyordu" yerine "agirliklama uniform
olarak daha kotu performans veriyor ve en kucuk kosulda sinyali kaybediyor".

## 3. Eslesmis 20-bolge cekirdegi (agirliksiz)

| yuk (ng) | n | Spearman | dMAE | P |
|---|---|---|---|---|
| 250  | 20 | 0.164 | 0.054 | 0.0796 |
| 500  | 20 | -0.224 | 0.027 | 0.1940 |
| 1000 | 20 | 0.134 | 0.016 | 0.2637 |
| 2000 | 20 | 0.405 | 0.034 | 0.1294 |

**0/4.** Cekirdek 20 bolge zaten tam oznitelikliydi, bu yuzden 44-bolge seti
bu sonucu DEGISTIRMEDI. Rastgele 20-altkume referansinda cekirdek %3
percentile'da (agirliksiz).

## 4. Sentez — makalenin nihai iddiasi

- Yerel yapisal tanimlayicilar, dort substrat yukunun DORDUNDE de bolgeye
  ozgu CMLase yanitina dair permutasyon-kontrollu ongorucu bilgi tasiyor.
- Iliski bolgeler arasinda DUZGUN DAGILMIYOR: dort kosulda da guvenle
  olculen 20 lizinde sinyal hicbir kosulda yok ve bu altkume rastgele
  dagilimin %3'unde. Nedeni belirlenemiyor; olcum guveni yaniti ongormuyor
  ve yanit dagilimlari ayirt edilemez.
- Ters-varyans agirliklama uniform olarak daha kotu ve etkin orneklemi
  nominalin ~%30'una indiriyor; oran-temelli PTM yanit etiketlerinde
  varsayilan olarak kullanilmamali.
