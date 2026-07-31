# Stage 4 v0.5.5 — AlphaFold Canonical/Fragment Selection Fix

## Düzeltilen sorunlar

1. AlphaFold DB API aynı UniProt sorgusu için kanonik, isoform ve farklı veri
   kümelerine ait birden fazla model döndürebilir.
2. Eski seçim algoritması hedef siteyi içeren en kısa modeli seçtiği için
   kanonik dizideki siteye bir isoform yapısı atanabiliyordu.
3. 25 Haziran 2026 sonrasında yeni API alanları `sequenceStart`,
   `sequenceEnd`, `sequence` ve `modelEntityId` kullanılmaktadır.
4. Uzun protein fragmentlerinde koordinat numarası ile UniProt numarası
   çakışırsa doğrudan numaraya öncelik verilmesi yanlış kalıntıyı seçebilirdi.

## Yeni güvenlik katmanları

- Sorgulanan kanonik accession ile AlphaFold model accession'ı eşleştirilir.
- API model dizisi kanonik UniProt dizisiyle ve hedef kalıntıyla karşılaştırılır.
- Isoform model, kanonik model varken seçilmez.
- `sequenceStart/sequenceEnd` ve eski alanlar birlikte desteklenir.
- Fragmentlerde kalıntı önce API dizisi içindeki yerel indeksle eşlenir.
- Çıktıya `alphafold_model_accession`, `alphafold_sequence_match` ve
  `alphafold_site_residue_match` denetim sütunları eklenir.

## Uygulama

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate

unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Stage4_AlphaFold_Fix_v0.5.5.zip -d .
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

## Eski yapısal cache'i temizleme

Metadata ve indirilen mmCIF dosyaları korunabilir; yalnız eski model seçim
sonuçları silinmelidir.

```bash
rm -rf data/interim/stage4/alphafold_cache/features_by_protein
rm -f data/processed/stage4/stage4_candidate_structure_features.csv
rm -f data/processed/stage4/stage4_scored_candidates.csv
rm -f results/stage4/stage4_excluded_structure.csv
```

## Yapıları yeniden üretme

```bash
set -o pipefail

python -u scripts/run_stage4.py structures   --workers 2   2>&1 | tee results/stage4_structures_v055_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

`--overwrite-structures` gerekli değildir. Mevcut metadata ve doğru mmCIF
dosyaları yeniden kullanılır; eksik kanonik modeller otomatik indirilir.

## Denetim

```bash
python - <<'PY'
import pandas as pd
import re

f = pd.read_csv(
    "data/processed/stage4/stage4_candidate_structure_features.csv"
)
x = pd.read_csv(
    "results/stage4/stage4_excluded_structure.csv"
)

all_rows = pd.concat([f, x], ignore_index=True, sort=False)
isoform = all_rows["alphafold_entry_id"].fillna("").str.match(
    r"^AF-[A-Z0-9]{6,10}-\d+-F\d+"
)

print("Yapı-uygun:", len(f))
print("Yapı-dışlanan:", len(x))
print("Seçilen isoform model:", int(isoform.sum()))
print(x["structure_quality_reason"].value_counts(dropna=False))

assert int(isoform.sum()) == 0
assert not x["structure_quality_reason"].fillna("").str.startswith(
    "structure_residue_is_"
).any()

print("\nPASS: canonical AlphaFold models are used")
PY
```
