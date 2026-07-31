# Stage 4 v0.5.7 — Full Atlas Structure Completion Fix

## Düzeltilen sorunlar

1. Quick koşuda bir protein için yalnız seçilen bölgeler cache'e yazıldığında, full koşudaki aynı proteine ait ek bölgeler işlenmiyordu.
2. Geçici DNS/bağlantı hataları `structure_processing_error:*` olarak cache'e yazılıyor ve sonraki koşularda yeniden denenmiyordu.

v0.5.7 cache'i candidate ID düzeyinde denetler. Başarılı satırları korur; eksik veya geçici hata içeren satırları otomatik yeniden işler.

## Uygulama

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Stage4_Structure_Completion_Fix_v0.5.7.zip -d .
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

## Yeniden çalışma

`features_by_protein` cache klasörünü silmeyin. Yama yalnız eksik/yeniden denenebilir bölgeleri bulur.

```bash
python -u scripts/run_stage4.py structures --workers 2 \
  2>&1 | tee results/stage4_structures_full_v057_console.txt
```

Aşağıdaki iki sayaç sıfır olana kadar aynı komut güvenle yeniden çalıştırılabilir:

```bash
python - <<'PY'
import pandas as pd
x = pd.read_csv('results/stage4/stage4_excluded_structure.csv')
r = x['structure_quality_reason'].fillna('').astype(str)
print('Retryable network/processing:', r.str.startswith('structure_processing_error:').sum())
print('Cache-incomplete rows:', r.eq('structure_not_eligible').sum())
PY
```

İki değer de 0 olduktan sonra `score --refit-model`, `prioritize` ve `audit` yeniden çalıştırılmalıdır.
