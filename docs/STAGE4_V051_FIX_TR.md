# Stage 4 v0.5.1 Excel Parser Düzeltmesi

## Sorun

v0.5.0, Supplementary Excel dosyalarında gerçek sütun başlıklarının ilk satırda
olduğunu varsayıyordu. Nature ek dosyalarında açıklama satırları veya çok katmanlı
başlıklar bulunduğunda accession ve site sütunları bulunamıyordu.

## Düzeltme

Yeni parser ilk 40 satırı tarar, gerçek başlık satırını bulur, gerekirse iki başlık
satırını birleştirir ve sütun adları yetersizse içerikten UniProt accession/site
sütunlarını çıkarır. Parse başarısız olsa bile diagnostics CSV yazılır.

## Yeniden çalıştırma

```bash
python scripts/run_stage4.py prepare --mode quick
cat results/stage4/stage4_parse_diagnostics.csv
```

Prepare başarılıysa:

```bash
python scripts/run_stage4.py annotate
python scripts/run_stage4.py structures --workers 2
python scripts/run_stage4.py score
python scripts/run_stage4.py prioritize --top-n 50
python scripts/run_stage4.py audit
```
