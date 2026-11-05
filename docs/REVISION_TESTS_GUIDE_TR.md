# CML-RepairNet Hakem Revizyon Testleri — v0.8.0

Bu yama mevcut Stage 6 v0.7.2 sonuçlarını **silmez veya üzerine yazmaz**. Yeni çıktılar
`results/revision_tests/` altında oluşturulur.

## 1. Proje kökü

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
```

## 2. Önce mevcut durumu yedekleme

```bash
chmod +x scripts/backup_before_revision.sh
bash scripts/backup_before_revision.sh
```

Standart yedek; kaynak kodu, testleri, yapılandırmaları, belgeleri, `results/`, `models/`
ve `data/processed/` klasörünü içerir. Büyük ve yeniden indirilebilir ara önbellekleri de
saklamak için:

```bash
bash scripts/backup_before_revision.sh --include-caches
```

Başarılı işlem sonunda `backups/CML_RepairNet_pre_revision_*.tar.gz` ve SHA-256 dosyası
oluşur. `PASS` görülmeden revizyon analizlerine geçmeyin.

## 3. Yamayı kurduktan sonra regresyon testleri

```bash
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

Mevcut 69 teste bu yamanın 3 testi eklenir; mevcut yerel pakete göre toplamın en az 72
olması beklenir.

## 4. Girdi denetimi

```bash
python scripts/run_revision_tests.py doctor
```

Kontrol edilenler:

- 33 BSA etiketinin ve üç girdi CSV'sinin bire bir eşleşmesi,
- hedef değişkenin varlığı,
- transfer-safe ailelerden B-factor sütunlarının gerçekten çıkarılması,
- dizi, yapı ve proxy ailelerindeki özellik sayıları.

## 5. T1 — hızlı özellik ailesi ablasyonu

```bash
python -u scripts/run_revision_tests.py ablation \
  --mode quick \
  --models extra_trees ridge lightgbm \
  --groups 4 \
  --n-jobs 4 \
  2>&1 | tee results/revision_tests/t1_quick_console.txt
```

Quick koşu, Ward anchor + 5 benzersiz random-init uzamsal partisyon kullanır.

Ana aileler:

1. `sequence`
2. `structure_full`
3. `structure_transfer_safe`
4. `combined_full`
5. `combined_transfer_safe`
6. `combined_transfer_safe_plus_proxy`
7. `null_mean`

Pozitif `mean_delta_mae`, aday ailenin daha iyi olduğunu gösterir.

## 6. T2 — permuted-structure negatif kontrolü

Önce 200 iterasyonluk hızlı kontrol:

```bash
python -u scripts/run_revision_tests.py permutation \
  --iterations 200 \
  --feature-set combined_transfer_safe \
  --model extra_trees \
  --groups 4 \
  --partition-method ward \
  --n-jobs 4 \
  2>&1 | tee results/revision_tests/t2_quick_console.txt
```

Her outer fold'da yalnız **eğitim satırlarındaki** yapı özellikleri, ortak bir satır
permütasyonuyla karıştırılır. Böylece yapı sütunlarının iç korelasyonu korunur; site-yapı
eşleşmesi bozulur. Test fold'una ve dizi özelliklerine dokunulmaz.

## 7. Tam analiz

Hızlı sonuçlar ve dosya denetimleri sorunsuzsa:

```bash
python -u scripts/run_revision_tests.py all \
  --mode full \
  --partitions 50 \
  --groups 3 4 5 \
  --models extra_trees ridge lightgbm \
  --permutations 1000 \
  --permutation-model extra_trees \
  --n-jobs 4 \
  2>&1 | tee results/revision_tests/revision_full_console.txt
```

Bu tam çağrı 3/4/5 uzamsal grup için Ward anchor ve her grup sayısında 50 benzersiz
partisyon üretir. Hesap yükü yüksekse önce yalnız `--groups 4` ile tamamlayın.

## 8. Ana çıktılar

```text
results/revision_tests/
├── revision_doctor.json
├── revision_t1_ablation_partition_metrics.csv
├── revision_t1_ablation_oof_predictions.csv
├── revision_t1_ablation_summary.csv
├── revision_t1_pairwise_deltas.csv
├── revision_t1_pairwise_summary.csv
├── revision_t1_partition_manifest.csv
├── revision_t1_feature_manifest.csv
├── revision_t2_structure_permutation_iterations.csv
├── revision_t2_structure_permutation_summary.csv
├── revision_t2_real_structure_oof_predictions.csv
└── revision_tests_report.md
```

## 9. Önceden belirlenmiş karar kuralları

- `sequence → combined_transfer_safe` delta-MAE çoğu partisyon için pozitif olmalı.
- Contact-based proxy yalnız `combined_transfer_safe` modelini istikrarlı biçimde iyileştirirse
  korunmalı.
- B-factor katkısı yalnız BSA içi mekanistik açıklama olarak raporlanmalı; AlphaFold aktarımına
  taşınmamalı.
- Permütasyon p-değeri destekleyici değilse “structure-aware predictive signal” iddiası
  yumuşatılmalı.
- Hiçbir sonuç BSA dışı deneysel kalibrasyon olarak sunulmamalı.
