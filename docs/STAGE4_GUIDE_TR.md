# CML-RepairNet Stage 4 — İnsan/Fare CML Aday Atlası

## Amaç

Stage 4, yayımlanmış CML proteomik bölgelerini insan ve fare proteinlerine ait
AlphaFold yapılarıyla eşleştirir ve her bölge için:

- transfer-safe CML onarılabilirlik tahmini,
- model konsensüsü,
- conformal tahmin aralığı,
- yapısal OOD yüzdeliği,
- AlphaFold yapı güveni,
- deneysel kaynak kanıtı,
- Pareto cephesi ve portföy sırası

üretir.

Çıktılar **hesaplamalı önceliklendirmedir**; insan/fare proteinlerinde deneysel
onarım doğrulaması anlamına gelmez.

## Kritik transfer güvenliği

BSA modelindeki kristalografik B-factor ile AlphaFold pLDDT aynı fiziksel büyüklük
değildir. Stage 4 bu iki değişkeni birbirine dönüştürmez ve aynı model girdisi gibi
kullanmaz.

Atlas modelinden şu iki B-factor özelliği çıkarılmıştır:

- `bfactor_ca`
- `neighbor_mean_bfactor_8a`

pLDDT yalnızca yapı kalite filtresi ve applicability-domain değerlendirmesinde
kullanılır.

## Veri kaynakları

- Di Sanzo et al. (2021), Supplementary Data 1 ve 2
- UniProt REST API
- AlphaFold Protein Structure Database API
- Stage 1 BSA etiket ve yapı özellikleri

İndirme URL'leri `src/cml_repairnet/atlas/sources.py` içinde sürüm kontrollüdür.
Her dosya için SHA-256 manifesti yazılır.

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Stage4_Patch_v0.5.0.zip -d .
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

## Ön kontrol

```bash
python scripts/run_stage4.py doctor
```

Gerekli dosyalar:

```text
data/processed/bsa_labels.csv
data/processed/bsa_structure_features.csv
```

## Quick çalışma

Quick mod, türler arasında dengeli yaklaşık 25 protein/60 bölge seçer. Amaç API,
struktur çıkarımı, model skorlama ve Pareto hattının birlikte çalıştığını doğrulamaktır.

```bash
set -o pipefail
python -u scripts/run_stage4.py all \
  --mode quick \
  --workers 2 \
  2>&1 | tee results/stage4_quick_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

## Parça parça çalıştırma

```bash
python scripts/run_stage4.py download
python scripts/run_stage4.py prepare --mode quick
python scripts/run_stage4.py annotate
python scripts/run_stage4.py structures --workers 2
python scripts/run_stage4.py score
python scripts/run_stage4.py prioritize --top-n 50
python scripts/run_stage4.py audit
```

Springer Nature indirmesi kurumsal ağ/proxy nedeniyle başarısız olursa, dosyaları
elle aşağıdaki adlarla yerleştirin:

```text
data/raw/stage4/cellular_cml_sites.xlsx
data/raw/stage4/mouse_organ_cml_sites.xlsx
```

Supplementary tablo sütunları otomatik tanınmazsa açık CSV şablonu kullanılabilir:

```bash
python scripts/run_stage4.py all \
  --mode quick \
  --candidate-csv examples/stage4_candidate_template.csv
```

## Full çalışma

Quick audit PASS olduktan sonra:

```bash
set -o pipefail
python -u scripts/run_stage4.py all \
  --mode full \
  --workers 4 \
  --top-n 100 \
  2>&1 | tee results/stage4_full_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

İşlem kaldığı yerden devam eder:

- UniProt kayıtları JSON cache olarak tutulur.
- AlphaFold metadata ve mmCIF dosyaları cache edilir.
- Her protein için özellik CSV'si ayrı yazılır.
- `--overwrite-structures` verilmedikçe tamamlanmış proteinler yeniden işlenmez.

## Ana çıktılar

```text
models/stage4/stage4_transfer_safe_model_bundle.joblib
models/stage4/stage4_transfer_safe_model_card.json
models/stage4/stage4_transfer_safe_oof_predictions.csv

data/processed/stage4/stage4_candidate_structure_features.csv
data/processed/stage4/stage4_scored_candidates.csv

results/stage4/atlas/cml_repairability_atlas_all_sites.csv
results/stage4/atlas/cml_repairability_atlas_in_domain.csv
results/stage4/atlas/cml_repairability_pareto_front.csv
results/stage4/atlas/cml_repairability_priority_portfolio.csv

results/stage4/stage4_excluded_sequence.csv
results/stage4/stage4_excluded_structure.csv
results/stage4/stage4_audit.json
results/stage4/stage4_audit.md
```

## Model portföyü

### Ana transfer tahmini

Nested-pruned Extra Trees, B-factor özellikleri çıkarılarak tüm 33 BSA bölgesi
üzerinde yeniden eğitilir.

### İkincil sıralama

Clipped ölçüm-belirsizliği ağırlıklı LightGBM.

### Yorum modeli

Ağırlıksız Ridge.

### Conformal belirsizlik

Transfer-safe Extra Trees'in 4-fold Ward spatial OOF artıklarından hesaplanır.
Tahmin aralıkları `[0,1]` sınırlarına kırpılır.

### OOD

Transfer-safe BSA özellik uzayında standartlaştırılmış KNN mesafesi kullanılır:

- `<0.95`: `in_domain`
- `0.95–0.99`: `borderline`
- `≥0.99`: `out_of_domain`

## Pareto amaçları

Maksimize edilenler:

- Extra Trees repairability
- Model konsensüs sırası
- Kaynak kanıt skoru
- Yapı kalite skoru

Minimize edilenler:

- OOD yüzdeliği
- Tahmin aralığı genişliği
- Modeller arası anlaşmazlık

Pareto front, ağırlıklı skorun alternatifi değil; çok amaçlı ana önceliklendirme
çıktısıdır. `weighted_priority_score` yalnız pratik sıralama kolaylığı sağlar.

## Yayın güvenlik dili

Kullanılabilir:

> The atlas provides computational prioritization of experimentally observed CML
> sites under explicit structural domain-shift and uncertainty controls.

Kullanılmamalı:

> The atlas validates repairability across the human or mouse proteome.
