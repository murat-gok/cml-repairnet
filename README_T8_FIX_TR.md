# CML-RepairNet Revision Phase II v0.11.1 — T8 düzeltmesi

Bu küçük yama, başarılı g:Profiler sonuçlarından `term_robustness` tablosu üretilirken oluşan
`KeyError: ['significant_background_list'] not found in axis` hatasını düzeltir.

## Değişiklik

`significant_background_list` sütunu ilk çevrimiçi çalıştırmada henüz mevcut olmadığından,
sütun kaldırma işlemi `errors="ignore"` ile idempotent hale getirildi. Bilimsel analiz,
arka planlar, eşikler ve sonuç hesapları değiştirilmedi.

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Phase2_T8_Fix_v0.11.1.zip -d .
python -m pytest -q tests/test_revision_phase2.py
```

## T8'i tekrar çalıştırma

Mevcut g:Profiler cache dosyaları güvenle yeniden kullanılabilir:

```bash
set -o pipefail
python -u scripts/run_revision_phase2.py enrichment --controls-per-query 5 \
  2>&1 | tee results/revision_phase2_t8_console.txt
echo "EXIT_CODE=${PIPESTATUS[0]}"
```
