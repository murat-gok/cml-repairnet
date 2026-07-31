# CML-RepairNet Stage 5 v0.6.0 — Biyolojik Yorumlama Rehberi

## Amaç

Stage 5, Stage 4'te dondurulan insan/fare CML onarılabilirlik atlasını yeni bir
makine öğrenmesi yarışmasına sokmaz. Amaç, dondurulmuş tahminleri biyolojik olarak
yorumlamak ve makaleye dönüştürülebilir çıktılar üretmektir.

Analizler:

1. UniProt protein ailesi, işlev, yolak, GO, InterPro ve Pfam anotasyonları.
2. g:Profiler ile GO, Reactome, KEGG, WikiPathways, CORUM ve HPA zenginleştirmesi.
3. STRING v12.0 ile yüksek güvenli protein–protein etkileşim ağları.
4. İnsan–fare ortolog eşleme ve skor uyumu.
5. Yapısal/tahminsel/kanıt özelliklerine göre aday modülleri.
6. Makale tabloları, PNG/PDF şekiller ve otomatik yorumlama raporu.

## Temel bilimsel güvenlik kararı

Zenginleştirme analizlerinde tüm genom arka planı kullanılmaz. Her tür için
Stage 4'te **in-domain** kabul edilen genler özel arka plan olarak kullanılır.
Böylece CML adaylarının seçilme sürecinden doğan aday-evreni yanlılığı azaltılır.

Stage 5 sonuçları biyolojik hipotez ve işlevsel bağlam üretir. Stage 4 tahminlerinin
bağımsız deneysel doğrulaması değildir.

## Kurulum

Yamayı proje köküne açtıktan sonra:

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

Stage 4 v0.5.7 sonrasında 46 test bulunduğu için beklenen toplam:

```text
55 passed
```

## Ön kontrol

```bash
python scripts/run_stage5.py doctor
```

Beklenen:

```json
{
  "portfolio": {"exists": true},
  "pareto": {"exists": true},
  "atlas_in_domain": {"exists": true},
  "ready": true
}
```

## Önerilen kontrollü çalıştırma

### 1. UniProt anotasyonları

```bash
python -u scripts/run_stage5.py annotate \
  2>&1 | tee results/stage5_annotate_console.txt
```

### 2. İşlevsel zenginleştirme

```bash
python -u scripts/run_stage5.py enrich \
  2>&1 | tee results/stage5_enrich_console.txt
```

Her tür için dört sorgu çalışır:

- `portfolio`
- `pareto`
- `evidence_rich`
- `high_repairability`

İnsan ve fare analizleri ayrı yürütülür; türler aynı zenginleştirme testinde
karıştırılmaz.

### 3. STRING ağları

```bash
python -u scripts/run_stage5.py network \
  2>&1 | tee results/stage5_network_console.txt
```

Varsayılanlar:

- STRING v12.0 sabit API adresi
- functional network
- required score = 700
- eklenen dış düğüm = 0
- yalnız portföy ve Pareto proteinleri arasındaki ağ

### 4. İnsan–fare ortolog analizi

```bash
python -u scripts/run_stage5.py orthology \
  2>&1 | tee results/stage5_orthology_console.txt
```

### 5. Aday modülleri

```bash
python -u scripts/run_stage5.py modules \
  2>&1 | tee results/stage5_modules_console.txt
```

K sayısı 2–6 arasında silhouette skoru ile seçilir. KMeans sabit random seed ile
çalışır.

### 6. Şekiller ve rapor

```bash
python scripts/run_stage5.py figures
python scripts/run_stage5.py report
code results/stage5/stage5_biological_interpretation.md
```

## Tek komut

```bash
set -o pipefail
python -u scripts/run_stage5.py all \
  2>&1 | tee results/stage5_full_console.txt
echo "EXIT_CODE=${PIPESTATUS[0]}"
```

## Çevrimdışı yeniden üretim

İlk başarılı çalıştırmadan sonra API yanıtları şurada cache'lenir:

```text
data/interim/stage5/api_cache/
```

Aynı dış veritabanı anlık görüntülerini kullanarak:

```bash
python scripts/run_stage5.py all --offline
```

`--offline` modunda cache'i bulunmayan dış analizler açıkça
`offline_cache_missing` olarak raporlanır; aday modülleri, yerel tablolar ve rapor
yine üretilebilir.

## Ana çıktılar

### UniProt

```text
data/processed/stage5/stage5_uniprot_annotations.csv
data/processed/stage5/stage5_protein_family_summary.csv
```

### Zenginleştirme

```text
data/processed/stage5/stage5_gprofiler_enrichment_all.csv
data/processed/stage5/stage5_gprofiler_enrichment_leading.csv
results/stage5/stage5_enrichment_manifest.csv
results/stage5/stage5_enrichment_data_versions.csv
```

### STRING

```text
data/processed/stage5/stage5_string_edges.csv
data/processed/stage5/stage5_string_nodes.csv
data/processed/stage5/stage5_string_network_metrics.csv
data/processed/stage5/stage5_string_ppi_enrichment.csv
```

### Ortolog ve modüller

```text
data/processed/stage5/stage5_ortholog_mapping.csv
data/processed/stage5/stage5_cross_species_pairs.csv
data/processed/stage5/stage5_cross_species_summary.csv
data/processed/stage5/stage5_candidate_module_assignments.csv
data/processed/stage5/stage5_candidate_module_profiles.csv
```

### Makale çıktıları

```text
data/processed/stage5/stage5_manuscript_candidate_table.csv
results/stage5/figures/*.png
results/stage5/figures/*.pdf
results/stage5/stage5_biological_interpretation.md
results/stage5/stage5_audit.json
```

## Git notu

Geçici API cache'ini normal Git commit'ine eklememek için `.gitignore` dosyasına:

```text
data/interim/stage5/api_cache/
```

ekleyin. Nihai analizde kullanılan dış veritabanı sürümleri ve sorgu manifestleri
`results/stage5/` altında tutulur. Tam API cache'i gerekiyorsa Zenodo sürümüne ayrı
bir reproducibility archive olarak eklenebilir.
