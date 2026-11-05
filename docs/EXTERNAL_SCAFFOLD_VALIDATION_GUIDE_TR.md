# CML-RepairNet — BSA dışı protein-scaffold dış tutarlılık testi

## Bilimsel amaç

Bu analiz, BSA'da eğitilmiş dondurulmuş transfer-safe modelin, Trabosh ve arkadaşlarının
kazein, hemoglobin ve kollajen için bildirdiği **bulk CML ELISA** sıralamasını yeniden üretip
üretemediğini sınar.

Bu çalışma:

- modeli yeniden ayarlamaz,
- dış protein sonuçlarını hiperparametre seçiminde kullanmaz,
- BSA dışı verileri eğitim setine eklemez,
- site-resolved external validation iddiasında bulunmaz.

Kaynak Supplementary Table 4'te sıralama her üç modifikasyon düzeyinde de aynıdır:

```text
casein > hemoglobin > collagen
```

## Temsil paneli

Kaynak makale ticari proteinlerin kesin tür/izoform ve katalog bilgilerini vermediği için,
önceden tanımlı bovine canonical panel kullanılır:

- Casein: P02662, P02663, P02666, P02668
- Hemoglobin: P01966, P02070
- Type I collagen: P02453, P02465

Casein için hem eşit-bileşen hem de tipik whole-casein bileşim ağırlıkları; hemoglobin ve
collagen için stoichiometric ağırlıklar raporlanır. Bu panel doğrulanmış reagent kimliği değildir.

## Birincil analiz

```text
model                 = frozen Extra Trees
quality               = site and neighborhood pLDDT >= 70
applicability domain  = OOD percentile < 0.95
component weighting   = stoichiometric
site aggregation      = RSA-weighted mean
endpoint              = exact recovery of casein > hemoglobin > collagen
```

Ek olarak relaxed kalite, all-domain, eşit bileşen, ortalama, component medianı ve en
erişilebilir dörtte birlik dilim duyarlılıkları raporlanır.

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate

unzip -o \
  /mnt/c/Users/user/Downloads/CML_RepairNet_External_Scaffold_Validation_v0.10.0.zip \
  -d .

python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

## Çalıştırma

```bash
mkdir -p results/external_scaffold_validation
set -o pipefail

python -u scripts/run_external_scaffold_validation.py all \
  --bootstrap-iterations 2000 \
  2>&1 | tee results/external_scaffold_validation_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

İndirme sırasında geçici ağ hatası olursa aynı komut güvenle yeniden çalıştırılabilir; UniProt
ve AlphaFold dosyaları cache altında korunur.

## Ana çıktılar

```text
results/external_scaffold_validation/
├── external_scaffold_doctor.json
├── external_scaffold_component_manifest.csv
├── external_scaffold_structure_features.csv
├── external_scaffold_site_predictions.csv
├── external_scaffold_aggregates.csv
├── external_scaffold_rank_consistency.csv
├── external_scaffold_rank_bootstrap.csv
└── external_scaffold_report.md
```

## Sonuç yorumlama

Birincil test başarılı sayılabilmesi için:

1. Üç scaffold'un da strict/in-domain kapsamı bulunmalı,
2. sıralama casein > hemoglobin > collagen olmalı,
3. pairwise ordering accuracy 1.0 olmalı,
4. site-bootstrap exact-rank probability tercihen yüksek olmalı.

Strict kapsam yoksa sonuç başarısızlık değil, **modelin bu scaffold için applicability-domain
dışında olduğu** anlamına gelir. Relaxed sonuçlar yalnız tanısal olarak yorumlanır.
