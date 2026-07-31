# Stage 5 v0.6.1 Ortoloji ve Raporlama Düzeltmesi

Bu yama Stage 5 v0.6.0 sonuçlarında `ortholog_pairs=0` üreten g:Orth şema
ayrıştırma hatasını düzeltir.

## Kök neden

g:Orth yanıtında:

- `incoming`: giriş gen sembolü,
- `converted`: giriş/kaynak Ensembl gen kimliği,
- `ortholog_ensg`: hedef organizmanın Ensembl gen kimliği,
- `name`: hedef organizmanın gen sembolüdür.

v0.6.0, `converted` alanını hedef fare gen sembolü olarak kullandığı için insan
Ensembl kimliklerini fare portföyündeki gen sembolleriyle eşleştiremiyordu.

## Yama sonrası çalıştırma

```bash
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
python scripts/run_stage5.py --offline orthology
python scripts/run_stage5.py figures
python scripts/run_stage5.py report
```

Mevcut g:Orth JSON cache'i kullanılacağı için dış API çağrısı gerekmez. Cache yoksa
`--offline` seçeneğini kaldırın.

## Beklenen kontroller

```bash
python - <<'PY'
import pandas as pd
m = pd.read_csv('data/processed/stage5/stage5_ortholog_mapping.csv')
p = pd.read_csv('data/processed/stage5/stage5_cross_species_pairs.csv')
s = pd.read_csv('data/processed/stage5/stage5_cross_species_summary.csv')
print(m[['source_gene','target_gene','source_ensembl','target_ensembl']].head())
print('Portfolio ortholog pairs:', len(p))
print(s.to_string(index=False))
assert m['target_gene'].notna().sum() > 0
assert not m['target_gene'].astype(str).str.startswith('ENSG').any()
assert len(p) > 0
PY
```

Ayrıca zenginleştirme figüründe tür, sorgu ve veri kaynağı etiketleri gösterilir.
Modül silhouette skoru audit dosyasına eklenir ve düşük ayrışma halinde modüller
`exploratory_low_separation` olarak işaretlenir.
