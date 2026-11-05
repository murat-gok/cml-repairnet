# CML-RepairNet Hakem Revizyonu — T1b/T2-full/T3/T4

Bu yama önceki `results/revision_tests/` çıktısını değiştirmez. Yeni sonuçlar varsayılan olarak
`results/revision_confirmatory/` klasörüne yazılır. Çalışmalar tamamen in silico'dur.

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate

unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Reviewer_Confirmatory_Tests_v0.9.0.zip -d .
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

## Önerilen sıra

### 1. T1b — yapı, B-factor ve esneklik proxy ayrıştırması

```bash
set -o pipefail
python -u scripts/run_revision_confirmatory.py ablation \
  --groups 4 \
  --random-partitions 5 \
  --models extra_trees lightgbm ridge \
  --n-jobs 4 \
  2>&1 | tee results/revision_confirmatory_t1b_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

Yeni doğrudan karşılaştırmalar:

- `structure_transfer_safe` → `structure_transfer_safe_plus_proxy`
- `structure_transfer_safe` → `structure_full`
- `structure_transfer_safe_plus_proxy` → `structure_full_plus_proxy`
- B-factor-only ve proxy-only diagnostik modeller

### 2. T3 — nested spatial model selection

```bash
set -o pipefail
python -u scripts/run_revision_confirmatory.py nested \
  --groups 4 \
  --random-partitions 5 \
  --models extra_trees lightgbm ridge \
  --fixed-model extra_trees \
  --n-jobs 4 \
  2>&1 | tee results/revision_confirmatory_t3_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

Her dış spatial fold için model yalnız dış-eğitim verisi içindeki inner GroupKFold ile seçilir.
Aynı partisyon üzerinde geriye dönük seçilen `naive_best` model de hesaplanarak seçim iyimserliği
nicelendirilir.

### 3. T4 — frozen model için paired cluster bootstrap

T1b tamamlandıktan sonra:

```bash
python -u scripts/run_revision_confirmatory.py bootstrap \
  --predictions results/revision_confirmatory/revision_t1_ablation_oof_predictions.csv \
  --output results/revision_confirmatory \
  --partition-id ward_g4 \
  --feature-set structure_transfer_safe \
  --model extra_trees \
  --iterations 20000 \
  --sign-flip-iterations 100000
```

Birincil güven aralığı spatial-group cluster bootstrap'tır. Sıradan site bootstrap yalnız duyarlılık
analizi olarak raporlanır.

### 4. T2-full — bütün altı partisyon üzerinde negatif kontrol

Önce 200 iterasyonlu bütün-partisyon denetimi:

```bash
set -o pipefail
python -u scripts/run_revision_confirmatory.py permutation \
  --groups 4 \
  --random-partitions 5 \
  --feature-set structure_transfer_safe \
  --model extra_trees \
  --iterations 200 \
  --checkpoint-every 25 \
  --n-jobs 4 \
  2>&1 | tee results/revision_confirmatory_t2_check_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

Denetim PASS olursa aynı komutu `--iterations 1000` ile yeniden çalıştırın. Komut mevcut checkpoint'i
okur ve yalnız eksik iterasyonları tamamlar. Baştan başlatmak için `--no-resume` kullanılır.

## Tek komutlu hızlı doğrulama

```bash
python -u scripts/run_revision_confirmatory.py all \
  --groups 4 \
  --random-partitions 5 \
  --models extra_trees lightgbm ridge \
  --permutations 200 \
  --bootstrap-iterations 2000 \
  --sign-flip-iterations 10000 \
  --n-jobs 4
```

Bu komut hızlı bütünlük denetimidir. Nihai tablolar için yukarıdaki ayrı ve daha yüksek iterasyonlu
komutlar kullanılmalıdır.

## Ana çıktılar

- `revision_t1_ablation_summary.csv`
- `revision_t1_pairwise_summary.csv`
- `revision_t2_full_partition_summary.csv`
- `revision_t2_full_aggregate_summary.csv`
- `revision_t3_nested_partition_metrics.csv`
- `revision_t3_nested_selection_frequency.csv`
- `revision_t3_nested_summary.csv`
- `revision_t4_primary_observed_metrics.csv`
- `revision_t4_bootstrap_summary.csv`
- `revision_t4_paired_tests.csv`
- `revision_confirmatory_report.md`
