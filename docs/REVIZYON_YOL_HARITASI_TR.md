# CML-RepairNet — Revizyon Yol Haritası (TR)

Durum: PEDS-26-0122 editoryal ön değerlendirmede reddedildi (Prof. R. E. Campbell).
Bu belge, bağımsız hakem değerlendirmesi + doğrulanmış literatür taraması + repo/kod
incelemesi sonrası çıkarılan uygulanabilir revizyon planıdır.

Tarih: 2026-11-02

---

## 0. Repo incelemesinden çıkan İKİ YENİ BULGU

Bunlar hakem raporunda yoktu; kodu ve veriyi okuyunca ortaya çıktı ve planın
en yüksek getirili iki maddesi bunlar.

### BULGU-1 (kritik fırsat): Bulk external veri 9 gözlem içeriyor, kod bunu 3'e indiriyor

`data/external/trabosh2026_bulk_protein_repair.csv` 3 scaffold x 3 modifikasyon
seviyesi = **9 gözlem** içeriyor (repair_fraction + SD + n=4 replike, Trabosh Supp. Table 4).

Gözlenen değerler:

| scaffold | low | mid | high |
|---|---|---|---|
| casein | 1.34 (SD 0.61) | 0.97 | 0.92 |
| hemoglobin | 0.89 | 0.87 | 0.82 |
| collagen | 0.80 | 0.65 | 0.60 |

Buradan iki sinyal okunuyor:
1. **casein > hemoglobin > collagen sıralaması HER ÜÇ seviyede de korunuyor** —
   yani tek bir sıralama değil, üç bağımsız tekrarlı sıralama var.
2. **Her scaffold içinde modifikasyon seviyesi arttıkça onarım oranı monoton
   DÜŞÜYOR** (3/3 scaffold). Bu, substrat yüklenmesi/doygunluk etkisine işaret
   eden, makalede hiç kullanılmamış bir bulgu.

Ancak `src/cml_repairnet/revision/external_scaffold.py::observed_summary`
(satır ~479) üç seviyeyi **scaffold başına median'a indiriyor** → gözlenen vektör
n=3'e düşüyor → `compare_rankings` n=3 üzerinde Spearman hesaplıyor.
**n=3'te Spearman istatistiksel olarak anlamsızdır** (yalnız 6 permütasyon, min p≈0.167).
Cross-scaffold testinin "çökmesinin" sebebi kısmen budur: test zaten güçsüz kurgulanmış.

> Aksiyon: A2. Bu tek başına makalenin external validation bölümünü
> "başarısız test" olmaktan çıkarıp gerçek bir bulguya dönüştürebilir.

### BULGU-2 (raporlama açığı): Ters-varyans ağırlıklama kodda var, makalede yok

`src/cml_repairnet/data/prepare_labels.py` (satır 80-84) her bölge için
`repair_variance_approx` ve `label_weight = (1/varyans) / median` üretiyor;
`models/stage25.py::build_sample_weights` bunu uyguluyor ve
`atlas/modeling.py`, `revision/confirmatory.py`, `revision/external_scaffold.py`,
`revision/feature_ablation.py` hepsi `model__sample_weight=` ile fit ediyor.

Ağırlıklar çok çarpık: `bsa_labels.csv`'de label_weight 0.33 ile 9.6 arasında
(~30x fark). Yani bazı bölgeler fit'e 30 kat daha fazla katkı veriyor.

**Manuscript Methods bölümünde bu ağırlıklamadan hiç söz edilmiyor.**
Bu bir bilimsel hata değil — tam tersi, doğru olan davranış — ama
açıklanmamış bir ağırlıklama şeması tekrarlanabilirlik açığıdır ve hakem
bunu "undisclosed modeling choice" olarak işaretler.

> Aksiyon: B1 + G1.

---

## A. Zorunlu bilimsel revizyonlar

### A1 — Konumlandırmayı yeniden yaz (EN KRİTİK, pazarlık yok)

**Sorun:** Kaynak makale (Trabosh ve ark. 2026) aynı 33 bölgede aynı yapısal
analizi zaten niteliksel olarak yapmış ve yayınlamış. Doğrulanmış alıntılar:
- "CMLase showing higher activity on sites with neighboring cysteines and reduced
  activity on sites followed by hydrophobic regions (Supplementary Fig. 14)"
- "CMLase may favor more hydrophilic, solvent accessible sites"
- "Activity weakly correlated with local protein flexibility (B-factor)
  (Supplementary Fig. 13)"
- Fig. 3C üzerinde aynı yapıyı (PDB 4F5S) kullanmışlar.

Makalenin 9 tanımlayıcısı (Cys sayısı 8 A, hidrofobik oran 8 A, aromatik sayı 8 A,
SASA/RSA, temas sayıları) **tam olarak bu üç hipotezin nicelleştirilmesi**.

**Yapılacak:**
- Introduction'da Trabosh Supp. Fig. 13-14'e AÇIKÇA atıf ver ve
  "kaynak çalışma bu ilişkileri niteliksel olarak gözlemledi; biz bunları
  nicelleştirip leakage-kontrollü olarak test ediyoruz" de.
- "First predictive framework" / "we develop a framework" iddialarını kaldır.
- Yeni çerçeve: *"nicel doğrulama + transfer sınırlarının haritalanması"*.
- Kaynak makaledeki iç çelişkiye değin: Fig. 3C açıklaması "no clear preference
  based on primary or secondary structure" derken Results metni belirli tercihler
  bildiriyor. Bizim nicel sonucumuz bu çelişkiyi çözüyor → bu gerçek bir katkı.

**Efor:** 1 gün (yalnız metin). **Zorunlu.**
**Etki:** Bunu yapmazsan her hakem yakalar ve makale ölür. Yaparsan savunulabilir
bir "confirmatory quantification" makalesi olur.

### A2 — External validation'i 9 gozlem uzerine yeniden kur  [TAMAMLANDI — SONUC NEGATIF]

**Durum: 2026-11-02 itibariyla uygulandi ve calistirildi.**
Modul: `src/cml_repairnet/revision/external_level_test.py` (additive; frozen
`external_scaffold` pipeline'ina DOKUNULMADI)
Script: `scripts/run_external_level_test.py`
Test: `tests/test_external_level_test.py` (21 iddia, hepsi gecti)
Cikti: `results/external_level_test/`

#### Onemli istatistiksel duzeltme

Ilk degerlendirmede iki ayri hipotez karistirilmisti. Ayrilmalari zorunlu:

- **H-OBS:** *Gozlenen* siralama seviyeler arasi tekrarlanabilir mi?
  Bagimsiz birim = seviyeler (3).
- **H-MOD:** *Yapi-temelli model* bu siralamayi geri getiriyor mu?
  Bagimsiz birim = scaffold'lar (3).

Model tek bir **seviyeden bagimsiz** siralama uretiyor. Dolayisiyla seviye
sayisini artirmak H-MOD'un gucunu ARTIRMAZ. 3 scaffold ile tam-siralama
testinde ulasilabilir en kucuk p degeri **1/3! = 0.167**'dir ve bu bir
bilgi-teorik tavandir; hicbir analiz bunun altina inemez.

Yani ilk mesajda verilen "p = 0.0046" degeri H-OBS'a aittir, H-MOD'a DEGIL.

#### Gercek sonuclar

**H-OBS — gozlenen taraf (olumlu):**
- 3/3 seviyede siralama ayni: casein > hemoglobin > collagen
- p(herhangi tutarli siralama) = 6x(1/6)^3 = 0.0278
- p(on-belirtilmis siralama) = (1/6)^3 = 0.0046
- UYARI: seviyeler ayni protein ve ayni enzim partisi; korelasyonlu
  pseudo-replike. Bu p degerleri kanitin UST SINIRI olarak raporlanmali,
  bagimsiz-replike cikarimi olarak degil.

**H-MOD — model tarafi (OLUMSUZ, kritik):**
- 192 konfigurasyonun 96'si 3 scaffold'u da iceriyor (digerleri eksik kapsam)
- Bu 96'nin yalnizca **7'si (%7.3)** gozlenen siralamayi veriyor
- Sans beklentisi **%16.7** → model siralamayi geri getirme yonunde
  **hicbir egilim gostermiyor**
- En sik ongoruler: hemoglobin>casein>collagen (39x),
  collagen>casein>hemoglobin (35x) — ikisi de yanlis
- **On-belirtilmis birincil konfigurasyon hic degerlendirilemiyor**:
  kollajende uygun lizin yok (`incomplete_scaffold_coverage`)
- Not: 96 konfigurasyon ayni pipeline'in ic ice varyantlari, yani yuksek
  korelasyonlu. %7.3 betimleyicidir, "sanstan anlamli olarak kotu" diye
  raporlanamaz. Eslesen 7 konfigurasyonu post-hoc secmek specification-search
  artefakti olur — YAPILMAMALI.

**Doz-yanit (yeni, betimleyici):**
- 3/3 scaffold'da modifikasyon seviyesi arttikca onarim monoton azaliyor
- Ortalama agirlikli egim: -0.0615 / seviye
- Isaret testi p = 0.125 → 3 scaffold ile taban bu; **anlamli degil**, sadece
  dusundurucu. Yapisal model seviyeden bagimsiz oldugu icin bu etkiyi
  ongormuyor.

**Veri kalitesi:**
- casein/low = 1.34 (SD 0.61) fiziksel olarak imkansiz (>%100 onarim).
  Medyan kullanimi bunu gizliyordu. Modul artik bu satiri isaretliyor,
  winsorize ediyor ve ters-varyans agirligi sayesinde fit'te baskin olmasini
  engelliyor. `--no-winsorize` ile duyarlilik kosusu mevcut.

#### Bunun makale icin anlami

Ilk cerceveleme ("A2 en zayif bolumu en guclu bolume cevirir") **yanlisti**.
Dogru okuma:

1. External transfer **gercekten basarisiz** ve bu artik nicel olarak
   belgelenmis durumda — "degerlendirilemedi" belirsizligi yerine
   "96 konfigurasyonun %7.3'u, sans %16.7" somut ifadesi var.
2. Bu, makaleyi bir **temsil-gecerliligi (representation validity) makalesi**
   yapar; A1 ve D2 ile birebir ortusur. Savunulabilir ve yayimlanabilir bir
   negatif bulgudur.
3. Guc tavani (p >= 1/6, 3 scaffold) acikca raporlanmali. Bu, "neden duzgun
   test etmediniz" sorusunu onceden kapatir ve dürüstlük puani kazandirir.
4. "Prediktif framework" iddiasi bu sonucla birlikte TAMAMEN birakilmali.
   Kalan savunulabilir iddia: yapi-yanit sinyali BSA ICINDE var (permutasyon
   kontrolu P=0.001), ve BSA disina transfer olmuyor.

**Etki:** Makalenin cercevesi degisiyor — "tahmin cercevesi" degil,
"BSA-ici yapi-yanit sinyali + nicelenmis transfer basarisizligi".
Bu daha kucuk ama savunulabilir bir makale.

### A3 — R-kare'yi özetten ve ana metinden çıkar

**Sorun:** n=33'te R-kare model değerlendirmesi için güvenilmez bir metrik
(PLoS ONE 2020, doi:10.1371/journal.pone.0229345 — açıkça uyarıyor).
Ayrıca 9 öznitelik / 33 gözlem = 3.7:1, literatürdeki 10:1 kuralının altında
(J Big Data 2025, doi:10.1186/s40537-025-01346-9).

**Yapılacak:** Özette ve Results'ta birincil metrikler = MAE + Spearman rho +
permütasyon P. R-kare yalnız supplementary tabloda, "not used for model
selection" notuyla.

**Efor:** 2 saat. **Zorunlu.**

### A4 — Başlığı ve kapsam iddiasını küçült

Mevcut: *"Structure-aware modeling of site-specific carboxymethyllysine repair:
a framework for prioritizing age-related protein damage"*

Sorunlar: "repair" fiilen ölçülen şey değil (ölçülen = CML sinyal azalması);
"framework for prioritizing age-related protein damage" taşınamayan bir iddia.

Öneri: *"Local structural context predicts site-specific CMLase substrate
response in bovine serum albumin: quantification and transfer limits"*
veya benzeri — "substrate response" + "transfer limits" içeren bir başlık.

**Efor:** 1 saat. **Zorunlu.**

---

## B. Metodolojik revizyonlar

### B1 — Ters-varyans ağırlıklamayı Methods'ta açıkla (BULGU-2)
- `prepare_labels.py` formülünü ver: delta-method varyansı
  `var_a/c^2 + a^2*var_c/c^4`, ardından `1/(var+floor)`, median'a normalize.
- Hangi weight_mode kullanıldığını belirt (`MODEL_WEIGHT_MODES`).
- Ağırlık dağılımını (min/median/max = 0.33/1.0/9.6) supplementary'ye koy.
- **Ağırlıksız duyarlılık analizi ekle:** ağırlıklar kapalıyken MAE ne oluyor?
  Hakem bunu mutlaka soracak.
**Efor:** 1 gün (kısa kod + metin). **Zorunlu.**

### B2 — Effective sample size'ı açıkça raporla
Ağırlıklar çarpık olduğu için gerçek ESS 33'ten düşük.
Kish ESS = (sum w)^2 / sum(w^2) hesapla ve raporla. Bu dürüstlük puanı kazandırır
ve "n=33 bile değil" eleştirisini sen önce dile getirmiş olursun.
**Efor:** 2 saat. **Zorunlu.**

### B3 — Nested selection sonucunu lehine çevir
Mevcut: nested MAE 0.203 > sabit model 0.189, penalty 0.013.
Bunu "model-selection uncertainty" diye gömmek yerine açık bir mesaj yap:
*"bu veri boyutunda daha fazla arama daha kötü genelleme veriyor"* —
bu, C-bölümündeki metaheuristik sorusunun da cevabı (bkz. Bölüm H).
**Efor:** 3 saat (metin). Opsiyonel ama tavsiye edilir.

---

## C. Ek analizler

### C1 — Trabosh'un kendi niteliksel iddialarını tek tek test et
Kaynak makale 3 somut iddia ediyor. Her birini ayrı ayrı nicel test et:
| Trabosh iddiası | Bizim testimiz |
|---|---|
| Komşu Cys → aktivite yüksek | `cys_count_8A` tek değişkenli etki + CI |
| Hacimli hidrofobik/aromatik kanat → düşük | `hydrophobic_fraction_8A`, `aromatic_count_8A` |
| B-factor zayıf korelasyon | zaten var (ablation) — Trabosh Supp. Fig. 13 ile karşılaştır |
| Dirençli bölgeler K131/K504/K573 | model bu 3 bölgeyi doğru düşük sıralıyor mu? |

Son satır özellikle değerli: Trabosh K131, K504, K573'ü <%5 onarım olarak
bildiriyor. Model bu üçünü en alt sıralara koyuyorsa bu **bağımsız, bölge düzeyli
bir doğrulamadır** ve şu an makalede yok.
**Efor:** 2 gün. **Zorunlu** (A1'in kanıt tabanı).

### C2 — Occurrence vs repairability ayrımını nicel göster
Literatürdeki tüm glikasyon tahmincileri (NetGlycate, GlyStruct, PredGly,
BERT-Kgly, iGly-IDN) *oluşumu* tahmin ediyor. Eğer bir occurrence skoru
üretilebiliyorsa (ör. GlyStruct benzeri öznitelikler), onun repair yanıtı ile
korelasyonunu göster. Düşük korelasyon = "bunlar farklı problemler" iddiasının
doğrudan kanıtı. Bu, makalenin EN savunulabilir novelty iddiası ve şu an metinde yok.
**Efor:** 2-3 gün. Tavsiye edilir, zorunlu değil.

### C3 — ESM embedding'leri karşılaştırma tabanı olarak kullan
Repoda `features/esm.py` ve `scripts/extract_esm_embeddings.py` var.
ESM tabanlı bir model, yapısal 9-öznitelikli modelden daha kötü performans
veriyorsa bu **yapısal bilginin diziden gelmediğini** gösterir — güçlü bir argüman.
**Efor:** 1-2 gün (kod zaten var). Tavsiye edilir.

---

## D. Biyolojik ve yapısal doğrulama

### D1 — Atlas'ı küçült, mekanizmayı büyüt
Mevcut: 1.268 bölge + Pareto + 100-bölge portföy + STRING + g:Profiler.
İnsan enrichment'ında korreksiyon sonrası HİÇBİR terim ayakta kalmıyor;
fare tarafında 3 geniş tema kalıyor. Yüksek bağlamlı 10 bölgenin 7'si fare
albümini — eğitim iskeletine yapısal olarak çok yakın, yani bağımsız doğrulama değil.
**Yapılacak:**
- İnsan g:Profiler ve STRING bölümlerini ana metinden çıkar → supplementary.
- Pareto/NSGA-II ve 100-bölge portföyü supplementary'ye indir.
- Ana metinde yalnız: 9 Tier-1 bölge + bunların yapısal gerekçesi + delivery sorunu.
- Fare albümini bölgelerinin **neden doğrulama sayılmadığını** açıkça yaz
  (şu an yazıyorsun — bunu koru, hatta öne çıkar).
**Efor:** 2 gün. **Zorunlu.**

### D2 — Kollajenin neden tamamen elendiğini mekanistik olarak işle
Kollajende hiç uygun lizin yok çünkü AlphaFold monomeri üçlü sarmalı
temsil etmiyor. Bunu bir başarısızlık olarak değil, **temsil geçerliliği
bulgusu** olarak yaz: "monomerik tahmin yapıları fibriler substratlar için
geçersiz bir öznitelik uzayı üretir". Bu, alan için gerçek bir uyarı.
**Efor:** 1 gün (metin). **Zorunlu.**

---

## E. Yeniden yapılandırma

### E1 — Introduction
- Trabosh'un niteliksel bulgularını ekle (A1).
- Eksik literatürü ekle: amadoriaz / fruktozil peptit oksidaz (FPOX) ve FN3K'nın
  **yapı temelli substrat özgüllüğü** çalışmaları. CMLase mekanistik olarak tam
  bu ailede (FAD-bağımlı oksidaz, glisin oksidaz iskeleti).
  Şu an bunlar hiç atıf almıyor → hakem için "literatürü bilmiyor" sinyali.
- DJ-1/PARK7 deglikaz tartışmasını **çekişmeli** olarak ekle
  (Richarme 2015 vs PARK7-glyoxalase karşı kanıtı). Bu, "AGE onarımı iddialarında
  temkinli ol" mesajını verir ve senin temkinli tonunu destekler.
- Occurrence vs repairability ayrımını kur (C2).

### E2 — Results
Sıralamayı değiştir: (1) yapı-yanıt sinyali + permütasyon kontrolü,
(2) Trabosh iddialarının nicel testi (C1), (3) sınır/compression davranışı,
(4) external 9-gözlem testi + doz-yanıt bulgusu (A2), (5) transfer sınırları ve atlas (küçültülmüş).

### E3 — Discussion
"Framework" dilini "quantification + transfer limits" diline çevir.
Negatif sonuçları (kollajen, enrichment, nested penalty) bir araya topla ve
"bu alanda neyin işe yaramadığı" başlığı altında güçlü bir bölüm yap.
**Efor:** 3-4 gün toplam. **Zorunlu.**

---

## F. Novelty güçlendirme — savunulabilir 3 iddia

Yalnız bu üçünü iddia et, fazlasını değil:
1. **Onarılabilirliğin (repairability) sürekli değişken olarak ilk denetimli
   regresyonu** — mevcut tüm literatür modifikasyon *oluşumunu* sınıflandırıyor.
   (Doğrulandı: NetGlycate 2006 → GlyStruct 2018 → PredGly 2019 → BERT-Kgly 2022
   → iGly-IDN 2024 zincirinin tamamı occurrence.)
2. **Trabosh/CMLase bölge düzeyli verisinin ilk ikincil analizi**, applicability
   domain ve dürüstçe raporlanmış negatif transfer sonucuyla.
   (Doğrulandı: literatürde başka reanaliz yok.)
3. **Transfer geçerliliği haritası:** monomerik yapı tanımlayıcılarının
   düzensiz/oligomerik/fibriler substratlara transfer olmadığının gösterilmesi.

Bunları zayıflatan ve KABUL EDİLMESİ gereken üç gerçek:
1. Kaynak makale aynı analizi niteliksel yapmış (A1).
2. Yöntem şablonu rutin — sistein ligandability (DeepCoSI 2022; PMC11040703'teki
   ağaç tabanlı modeller neredeyse birebir aynı öznitelikleri kullanıyor),
   lizin PTM tahmini (UbiSite, KA-predictor), NGOME deamidasyon.
3. n=33, tek protein; kendi cross-scaffold testi başarısız.

---

## G. Tekrarlanabilirlik

- **G1:** Ağırlıklama şemasını Methods'a yaz (B1).
- **G2:** Bibliyografik düzeltmeler:
  - Di Sanzo ve ark. için **yayınlanmış başlık** kullan:
    *"Mapping protein carboxymethylation sites provides insights into their role
    in proteostasis and cell proliferation"* (Nat Commun 2021;12:6743).
    Makalede şu an preprint/dataset başlığı da geçiyor.
  - Trabosh: "pages 5926" → **"article 5926"** (Nature Communications makale numarası kullanır).
- **G3:** `--seed` ve freeze manifest'leri zaten var; bunları Methods'ta
  tek cümleyle işaret et ve GitHub DOI'si (Zenodo) al.
- **G4:** `data/external/trabosh2026_bulk_protein_repair.csv`'nin Supp. Table 4'ten
  nasıl çıkarıldığını belgele (şu an source kolonu var ama dönüştürme adımı yok).

---

## H. Sezgisel (metaheuristik) optimizasyon temelli hibrit yöntem — DEĞERLENDİRME

Ayrı bölüm; sohbetteki cevabın kalıcı kaydı.

### Modelin kendisi için: HAYIR

1. **Kendi verin zaten tersini kanıtlıyor.** Nested adaptif model seçimi
   MAE'yi 0.189'dan 0.203'e YÜKSELTTİ (penalty 0.013, CI 0.006-0.020).
   Yani bu veri boyutunda arama uzayını genişletmek genellemeyi kötüleştiriyor.
   DE/GA/PSO ile öznitelik seçimi veya hiperparametre araması bu etkiyi
   büyütür, küçültmez — klasik winner's curse.
2. **n=33, 9 öznitelik.** Metaheuristik arama 2^9 öznitelik altkümesi x
   hiperparametre uzayında optimize ederken 33 noktayı defalarca kullanır.
   Literatür net: n<100'de regresyon zaten kırılgan, R-kare güvenilmez.
3. **Hakemlerin itirazı "model yeterince optimize değil" DEĞİL.** İtiraz
   (a) soruyu kaynak makale zaten cevaplamış, (b) n çok küçük, (c) transfer
   başarısız. Metaheuristik bunların hiçbirine dokunmuyor.
4. **Desk reject riskini artırır.** PEDS reddinin en olası nedeni scope
   uyumsuzluğu + n=33'lük çekirdek üstünde aşırı yöntem iskelesi
   ("methodological theater") algısı. Bir protein bilimi dergisinde
   33 nokta üzerinde DE-tuned ERT, sofistikasyon değil kırmızı bayrak okunur.
5. **Alan hizalaması aleyhine çalışır.** Tıp fakültesine geçiş hedefi
   açısından makaleyi optimizasyon makinesiyle yüklemek, onu bilgisayar
   mühendisliği tarafına geri çeker.

### Tek meşru kullanım: modeli değil, SONRAKİ DENEYİ optimize et — EVET

Burada gerçek ve savunulabilir bir kombinatoryal optimizasyon problemi var:

> *k protein/bölgeden oluşan bir deney paneli seç; öyle ki yapısal tanımlayıcı
> uzayının kapsanması, öngörülen yanıt çeşitliliği, tür dengesi ve
> hücresel erişilebilirlik birlikte maksimize edilsin; albümin-benzeri
> kolay vakalara çökmesin.*

Bu bir **design-of-experiments / maximin-D-optimal / submodüler kapsama**
problemi ve:
- Küçük n bir dezavantaj değil, çünkü burada model fit etmiyorsun, **karar** veriyorsun.
- NSGA-II zaten repoda (`atlas/prioritize.py`, `robustness/ranking.py`) ve
  `active_learning/select_candidates.py` iskeleti mevcut.
- Makalenin gerçek katkısıyla (hipotez önceliklendirme) birebir örtüşüyor.
- Literatürde CML onarımı için böyle bir panel tasarımı yok → gerçek novelty.

Bu, makaleyi *"prediktif framework"* (savunulamaz) yerine
*"CML onarım substrat paneli için deney tasarımı çerçevesi"* (savunulabilir)
olarak konumlandırır ve A1'deki küçültmeyi telafi eder.

**Özet: modeli optimize etme, sonraki deneyi optimize et.**

---

## I. Fazlı uygulama planı

### Faz 1 — Metin ve ucuz düzeltmeler (3-4 gün, kod yok veya çok az)
- A1 konumlandırma
- A3 R-kare
- A4 başlık
- G2 bibliyografik düzeltmeler
- B1/G1 ağırlıklama Methods metni
- B2 Kish ESS (2 saat kod)

Faz 1 sonunda makale "savunulabilir" hale gelir. Buradan sonrası performans artışı.

### Faz 2 — Kod gerektiren çekirdek güçlendirme (1.5-2 hafta)
- A2 external 9-gözlem yeniden kurgusu  ← EN YÜKSEK GETİRİ
- C1 Trabosh iddialarının nicel testi (K131/K504/K573 dahil)
- B1 ağırlıksız duyarlılık analizi
- D1 atlas küçültme
- D2 kollajen temsil geçerliliği bölümü

### Faz 3 — Opsiyonel güçlendirme (1-2 hafta)
- C2 occurrence vs repairability nicel ayrımı
- C3 ESM karşılaştırma tabanı
- H meşru optimizasyon: panel tasarımı modülü

### Faz 4 — Yeniden gönderim
- Dergi seçimi (ayrı çalışma: SCIE, APC'siz, >=2 Springer Nature)
- Cover letter'da PEDS reddine değinme; scope uyumunu ilk paragrafta kur

---

## J. Gerçekçi hedef beklentisi

- **Mevcut haliyle gönderim:** önerilmez.
- **Faz 1 sonrası:** Q4, sınırda Q3.
- **Faz 1 + Faz 2 sonrası:** Q3 gerçekçi; iyi scope uyumuyla Q2 sınırı.
- **Faz 3 dahil:** Q2 savunulabilir hale gelir, ama garanti yok.

Not: Site-çözünürlüklü external doğrulama verisi LİTERATÜRDE YOK
(Trabosh'un panel ölçümleri bulk CML ELISA). Dolayısıyla "gerçek external
site-level validation" bu makale için ulaşılamaz; A2'deki 9-gözlem kurgusu
bu kısıt altında yapılabilecek en güçlü testtir.

---

## K. BULGU-3 (2026-11-02, EN ONEMLI) — PXD070124'te kullanilmayan 3 kosul var

PRIDE PXD070124 incelendi (API + FTP). Ham veriyi yeniden islemeye GEREK YOK:
islenmis bolge duzeyli niceleme CSV olarak zaten depolanmis (8 kucuk dosya).

Depodaki dosyalar: `BSA_{0p25,0p5,1p0,2p0}_865_794_processed_filt_quant.csv`
Kolonlar: Protein, Residue #, Modification, Peptide Sequence, Confidence,
ve 3 aktif + 3 inaktif replike icin "% Abundance".

| kosul | toplam bolge | Confidence>90 |
|---|---|---|
| 0p25 | 35 | 22 |
| 0p5  | 41 | **33  <- makalede kullanilan** |
| 1p0  | 44 | 37 |
| 2p0  | 45 | 38 |
| **toplam gozlem** | | **130** |

`0p5` oldugu KESIN: bolge 65 icin aktif 9.98/8.63/9.43 ve inaktif
17.64/17.40/17.92 degerleri `data/processed/bsa_labels.csv` ile birebir ayni,
hesaplanan yanit 0.470 ~ repair_published 0.4615. Yani >90 guven filtresi
dogru uygulanmis ve tekrarlanabilir — ama **dort kosuldan yalnizca biri
kullanilmis.**

### Dogru okuma (gozlem sayisi != bagimsiz birim)

Dort kosul ayni proteinin farkli enzim:substrat oranlari, bagimsiz biyolojik
replike DEGIL. Dolayisiyla "n 33'ten 130'a cikiyor" demek A2'de yaptigim
hatanin tekrari olur. Gercek kazanclar:

1. **Bagimsiz birim (bolge) sayisi artar:** 33 -> birlesimde ~45.
2. **Bolge basina tekrarli olcum:** etiket kesinligi ciddi artar; ters-varyans
   agirliklama (B1) ancak o zaman tam anlamini kazanir.
3. **Merkezi bulgu 4 kez tekrarlanabilir hale gelir.** Makalenin ana iddiasi
   su an TEK kosula dayaniyor. Permutasyon-kontrollu yapi-yanit sinyali
   dort kosulda da ayakta kalirsa bu, elde edilebilecek en guclu ic
   dogrulamadir — ve makalede su an boyle bir sey YOK.
4. **Yeni ve ozgun bir soru acilir:** yerel yapi, bolgenin *konsantrasyon
   duyarliligini* ongoruyor mu? Bu farkli bir bagimli degisken, ayni verinin
   daha fazlasi degil. Literaturde yok.

Ornek: bolge 65 yaniti kosullara gore 0.327 / 0.470 / 0.203 / 0.200 —
belirgin konsantrasyon bagimliligi var (bulk duzeyinde A2'de bulunan
doz etkisiyle tutarli).

### Risk

Sinyal dort kosulda tekrarlanmazsa makale daha da zayiflar. Ama bunu
bilmek gerekiyor — hakem zaten sorabilir ("neden tek kosul?").

### Maliyet

Saatler, haftalar degil. 8 kucuk CSV indirilecek, mevcut feature pipeline
aynen kullanilacak. Ham MS yeniden isleme (MaxQuant/FragPipe) GEREKSIZ.

**Oneri: Faz 2'nin basina alinmali, A1'den bile once bakilmali — cunku
sonucu makalenin cercevesini belirler.**
