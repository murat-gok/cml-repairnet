# CML-RepairNet tamamen in silico yol haritası

> Dosya adı önceki sürümle uyumluluk için korunmuştur. Projede ıslak laboratuvar, hayvan, klinik
> çalışma veya yeni deneysel veri üretimi yoktur.

## 1. Kesin kapsam ve bilimsel iddia

CML-RepairNet'in hedefi, CML oluşumunu tahmin etmek değil; **CML oluşmuş bir lizinin CMLase ile
onarılabilirliğini** 0–1 aralığında tahmin etmektir. Model her tahmin için belirsizlik, OOD skoru ve
kanıt düzeyi üretir.

Proje şu iddiayı kurabilir:

> Yayımlanmış site-bazlı CMLase verileri, protein dizisi, yapı, dinamik ve protein-dil modeli
> temsilleriyle birleştirildiğinde, CML onarılabilirliğini açıklayan hesaplamalı bir çerçeve ve
> proteom ölçekli hipotez atlası üretilebilir.

Proje şu iddiayı kuramaz:

> İnsan proteomundaki tahmin edilen bölgelerin deneysel veya klinik olarak doğrulandığı.

## 2. Ana hipotezler

- **H1:** Yerel üç boyutlu erişilebilirlik ve paketlenme, yalnızca birincil diziden daha fazla sinyal
  taşır.
- **H2:** Tek bir statik yapı yetersizdir; konformasyon ansamblından üretilen özellikler tahmin
  kararlılığını artırır.
- **H3:** Dondurulmuş protein dil modeli temsili ile fiziksel yapı özelliklerinin birleşimi, tekil
  özellik ailelerini geçer.
- **H4:** Çok-fideliteli hesaplamalı skorlar, 33 etiketli bölgeden öğrenilen modele ek bilgi sağlar.
- **H5:** Belirsizlik, biyolojik önem, yapısal çeşitlilik ve hesaplama maliyetini birlikte optimize
  eden seçim, yalnızca tahmin skoruna göre sıralamadan daha dengeli bir sanal hedef portföyü üretir.

## 3. Aşama 0 — veri yönetişimi ve tekrar üretilebilirlik

### Çıktılar

- Veri manifestosu, URL/DOI/accession, sürüm ve SHA256 kayıtları
- Lisans ve yeniden dağıtım matrisi
- Conda ortamı, sabit rastgele tohumlar ve çalıştırma manifestosu
- Analiz protokolünün sonuçlar görülmeden önce dondurulmuş sürümü

### Go/no-go kriteri

33 site-bazlı BSA etiketi ile tekrar verileri güvenilir biçimde yeniden oluşturulamıyorsa ana model
çalışması başlatılmaz. Çelişkiler açık bir veri denetim raporunda çözülür veya sınır olarak yazılır.

## 4. Aşama 1 — BSA etiketlerinin yeniden yapılandırılması

1. Makale kaynak XLSX dosyasından aktif ve inaktif enzim tekrarlarını çıkar.
2. Yayın numaralarını UniProt P02769 konumlarına eşleştir ve her konumu `K` olarak doğrula.
3. Yayımlanmış onarım oranını, tekrar ortalamalarından hesaplanan alternatif etiketlerle karşılaştır.
4. Varyans, düşük kontrol bolluğu ve aykırı tekrar bayrakları üret.
5. PXD070124 işlenmiş ve mümkünse ham verilerinden bağımsız yeniden üretim yap.
6. Etiket sürümleri arasında duyarlılık analizi tanımla.

### Başarı ölçütleri

- 33/33 bölgenin tutarlı eşleşmesi
- Yayın sıralamasının yeniden üretimi
- Etiket farklarının normalizasyon, peptit birleştirme veya grafik hesaplama adımlarından hangisine
  bağlı olduğunun belgelenmesi

## 5. Aşama 2 — leakage-free bilimsel baz modeller

### Özellik aileleri

- ±5, ±10 ve ±15 amino asitlik dizi pencereleri
- Hidrofobiklik, yük, aromatiklik, hacim, esneklik ve motif göstergeleri
- SASA/RSA, residue depth, B-factor veya pLDDT
- 8/10/12 Å komşuluk ve paketlenme yoğunluğu
- İkincil yapı, yerel eğrilik ve yüzey yamaları
- Elektrostatik potansiyel ve lokal polarite özetleri

### Modeller

- Ortalama/medyan null model
- Ridge ve Elastic Net
- SVR
- Random Forest ve Extra Trees
- Gaussian Process Regression
- Bayesian Ridge veya heteroscedastic regresyon

### Değerlendirme

- Uzamsal GroupKFold
- Motif-kümesi dışlama
- Yapısal çevre-kümesi dışlama
- Nested hiperparametre ayarı
- Tüm ölçekleme, PCA, seçim ve imputasyon eğitim fold'u içinde

### Birincil metrik

MAE. İkincil metrikler: Spearman, RMSE, kalibrasyon hatası ve tahmin aralığı kapsaması.

## 6. Aşama 3 — protein dil modeli ve çok-modlu füzyon

- ESM-2 gibi önceden eğitilmiş modellerden hedef kalıntı ve yerel pencere gömmeleri çıkarılır.
- 33 örnek üzerinde uçtan uca fine-tuning yapılmaz.
- Boyut indirgeme yalnızca eğitim fold'u içinde uygulanır.
- Erken füzyon, geç füzyon ve stacking karşılaştırılır.
- Dizi-only, yapı-only, PLM-only ve birleşik modeller için ablation yapılır.
- Protein adı, mutlak konum veya veri kaynağı gibi kestirme yollar modele verilmez.

## 7. Aşama 4 — konformasyon ansamblı ve dinamik özellikler

Tamamen hesaplamalı olarak her etiketli bölge için yapı belirsizliği modellenir:

1. Deneysel BSA yapısı ile AlphaFold yapısı karşılaştırılır.
2. Eksik yan zincirler ve alternatif protonasyon durumları hazırlanır.
3. Normal mode analysis veya kısa MD ansamblı ile konformerler üretilir.
4. Her konformerde SASA, residue depth, cep açıklığı ve komşuluk hesaplanır.
5. Özelliklerin ortalaması, varyansı ve uç değerleri modele eklenir.

### Ana test

Tek yapı özellikleri ile ansambl özellikleri, aynı cross-validation bölmeleri üzerinde
karşılaştırılır. İyileşme yoksa pahalı simülasyon katmanı zorunlu tutulmaz.

## 8. Aşama 5 — çok-fideliteli hesaplamalı doğrulama

### Seviye 0: hızlı özellikler

Dizi, PLM, statik yapı ve basit erişilebilirlik.

### Seviye 1: ansambl erişilebilirliği

Normal mode/MD konformerleri, zaman-ortalama SASA, yerel esneklik ve sterik erişim.

### Seviye 2: CMLase–substrat geometrisi

Exact CrGO-897 dizisi ve güvenilir yapı modeli kamuya açık bilgilerden yeniden oluşturulabiliyorsa:

- CML-modifiye peptit/protein parçalarının ensemble docking'i
- Aktif merkeze giriş yolu, CML yönelimi ve sterik çakışma skorları
- Birden fazla docking motoru veya parametre setiyle consensus skor

Dizi güvenilir biçimde yeniden oluşturulamıyorsa bu katman ana analizden çıkarılır; varsayımsal enzim
modeliyle kesin sonuç üretilmez.

### Seviye 3: seçilmiş adaylarda enerji ve dinamik

- Enzim–peptit komplekslerinde kısa MD
- Temas sürekliliği, aktif-geometri korunumu ve MM/GBSA benzeri göreli skorlar
- Yalnızca mekanizma ve parametreler savunulabiliyorsa küçük QM veya QM/MM duyarlılık analizi

Bu skorlar deneysel gerçek olarak değil, **yüksek maliyetli hesaplamalı yardımcı kanıt** olarak
kullanılır.

## 9. Aşama 6 — yayımlanmış protein-düzeyi verilerle dış tutarlılık testi

Kaynak çalışmada CMLase aktivitesi BSA dışında kollajen, hemoglobin, kazein ve göz proteini özütü
gibi farklı protein iskeletlerinde protein-düzeyi ölçülmüştür. Uyumlu sayısal kaynak veriler elde
edilebilirse:

1. Her protein için CML'lenebilir/raporlanmış lizinler belirlenir.
2. Site tahminleri ağırlıklı protein-düzeyi skora dönüştürülür.
3. Tahmin edilen protein sıralaması yayımlanmış global CML azalmasıyla karşılaştırılır.
4. Farklı agregasyon kuralları için duyarlılık analizi yapılır.

Bu test site-bazlı dış doğrulama değildir; ancak tek protein üzerinde öğrenilen modelin farklı
protein iskeletleriyle **dış tutarlılığını** sınar.

## 10. Aşama 7 — insan ve fare CML aday atlası

Yayımlanmış CML proteomik çalışmalarındaki 1000'den fazla bölge onarım etiketi olarak değil,
**sanal tarama aday havuzu** olarak kullanılır.

Her aday için:

- Dizi ve ortolog bilgisi
- PDB veya AlphaFold yapısı ve güven puanları
- CML bölgesinin deneysel kanıt düzeyi
- Doku, hücre tipi ve subselüler konum
- Protein yarı ömrü ve ECM ilişkisi
- Yaşlanma/hastalık ve PTM çakışmaları
- CML-RepairNet ortalama tahmini
- Konformal tahmin aralığı
- OOD skoru ve veri-kapsam etiketi

Atlasın her satırı `benchmark`, `interpolation`, `near-OOD` veya `far-OOD` kanıt sınıfıyla sunulur.

## 11. Aşama 8 — çok amaçlı sanal hedef optimizasyonu

Amaç, laboratuvar deneyi seçmek değil; en değerli **hesaplamalı hedef portföyünü** oluşturmaktır.

Maksimize edilecekler:

- Tahmin edilen onarılabilirlik
- Biyolojik/yaşlanma önemi
- Model güveni
- Yapısal ve doku çeşitliliği
- Bağımsız veri desteği

Minimize edilecekler:

- OOD riski
- Yapı belirsizliği
- Adaylar arası fazlalık
- Yüksek-fideliteli hesaplama maliyeti

Greedy taban çizgisi; NSGA-II, MOEA/D ve uygun durumda Bayesian multi-objective optimization ile
karşılaştırılır. Pareto cephesi tek bir “en iyi” liste yerine alternatif hedef portföyleri verir.

## 12. Aşama 9 — in silico aktif öğrenme

Yeni deneysel etiket gelmeyeceği için klasik aktif öğrenme yapılmaz. Bunun yerine **multi-fidelity
active learning** uygulanır:

1. Düşük maliyetli model tüm adayları puanlar.
2. Belirsiz ve çeşitli adaylar daha pahalı ansambl/docking/MD katmanına seçilir.
3. Yüksek-fideliteli hesaplamalı skorlar yardımcı hedef olarak modele eklenir.
4. Maliyet–bilgi kazanımı eğrisi rastgele ve yalnızca belirsizlik seçimleriyle karşılaştırılır.

Bu yaklaşımın çıktısı daha iyi deney seçimi değil, sınırlı GPU/CPU bütçesinin daha verimli
kullanılmasıdır.

## 13. Aşama 10 — sağlamlık, negatif kontroller ve yanlışlanabilirlik

- Tam pipeline label-permutation testi
- Koordinat ve motif karıştırma kontrolleri
- Özellik ailesi çıkarma ablation'ları
- Etiket tanımı duyarlılığı
- AlphaFold–PDB yapı duyarlılığı
- Farklı random seed ve fold tanımlarında kararlılık
- Yüksek repairability ve düşük repairability sınıflarında hata analizi
- Benzerlik tabanlı nearest-neighbor null modeli
- Tahminlerin basit SASA veya pLDDT kısayoluna indirgenip indirgenmediği
- Calibration ve conformal coverage denetimi

## 14. Nihai yayın paketi

- Denetlenmiş 33-site benchmark veri seti
- Leakage-free model karşılaştırması
- Belirsizlik kalibrasyonu ve OOD analizi
- Çok-fideliteli hesaplamalı kanıt katmanı
- Yayımlanmış farklı protein iskeletleriyle dış tutarlılık analizi
- İnsan/fare CML Repairability Atlas
- Pareto-optimal sanal hedef portföyleri
- Açık kod, ortam dosyası, veri manifestosu ve yeniden üretilebilir raporlar

## 15. Önerilen makale başlığı

**CML-RepairNet: A Structure-Aware Multi-Fidelity Framework for Predicting Site-Specific Enzymatic
Repairability of Protein Carboxymethylation**

Daha ihtiyatlı alternatif:

**CML-RepairNet: Leakage-Controlled Computational Modeling and Proteome-Scale Prioritization of
CML Repairability**

## 16. 12 aylık hesaplamalı takvim

- **Ay 1:** Veri yeniden üretimi, etiket denetimi ve baz modeller
- **Ay 2:** Yapı özellikleri, leakage testleri ve ablation
- **Ay 3:** PLM gömmeleri ve çok-modlu füzyon
- **Ay 4:** Konformasyon ansamblı ve dinamik özellikler
- **Ay 5:** Protein-düzeyi dış tutarlılık analizi
- **Ay 6–7:** CML aday atlasının kurulması
- **Ay 8:** OOD ve konformal belirsizlik
- **Ay 9–10:** Çok-fideliteli sanal simülasyonlar
- **Ay 11:** Çok amaçlı optimizasyon ve nihai atlas
- **Ay 12:** Makale, kod sürümü, model kartı ve arşivleme
