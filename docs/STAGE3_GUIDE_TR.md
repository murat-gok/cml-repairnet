# CML-RepairNet Stage 3 Kullanım Rehberi

Stage 3 yeni bir model yarışması değildir. Stage 2.5 sonunda önceden sabitlenen
model portföyünün genellenebilirlik, mekanizma ve belirsizlik denetimidir:

- Birincil sürekli tahmin: ağırlıksız, nested-pruned Extra Trees
- Aday sıralama: clipped-weight LightGBM
- Doğrusal yorum: Ridge

## Analiz bileşenleri

1. Seeded K-means ile repeated spatial cross-validation
2. Stage 2 Ward bölmesinin anchor analizi
3. BSA Domain I/II/III leave-one-domain-out değerlendirmesi
4. Cross-fitted affine ve isotonic boundary calibration
5. Group-aware split-conformal tahmin aralıkları
6. KNN tabanlı OOD skoru
7. `neighbor_cysteine_count_8a` için nested drop-column analizi
8. Training-fold conditional model ile residual conditional permutation

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

## Quick koşu

```bash
set -o pipefail
python -u scripts/run_stage3.py \
  --mode quick \
  --conditional-repeats 200 \
  --n-jobs 1 \
  2>&1 | tee results/stage3_quick_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

Quick mod, 5 seeded K-means bölmesi ve bir Ward anchor kullanır.

## Full koşu

```bash
set -o pipefail
python -u scripts/run_stage3.py \
  --mode full \
  --conditional-repeats 1000 \
  --n-jobs 1 \
  2>&1 | tee results/stage3_full_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

Full mod, 50 farklı K-means seed'i ve bir Ward anchor kullanır; ana grup sayısı 4'tür.

## Exhaustive grup duyarlılığı

```bash
python -u scripts/run_stage3.py \
  --mode exhaustive \
  --conditional-repeats 1000 \
  --n-jobs 1 \
  2>&1 | tee results/stage3_exhaustive_console.txt
```

Bu mod 50 seed'i 3, 4 ve 5 spatial group için çalıştırır.

## Özelleştirme

```bash
python scripts/run_stage3.py --mode full --seeds 100 --groups 4
python scripts/run_stage3.py --mode full --groups 3 4 5
python scripts/run_stage3.py --target repair_from_means_clipped --mode full
```

## Ana çıktılar

```text
results/stage3/stage3_report.md
results/stage3/repeated_spatial/repeated_spatial_summary.csv
results/stage3/repeated_spatial/repeated_spatial_metrics.csv
results/stage3/domain_holdout/domain_metrics.csv
results/stage3/domain_holdout/domain_uncertainty.csv
results/stage3/domain_holdout/domain_uncertainty_summary.csv
results/stage3/cysteine_mechanism/cysteine_mechanism_summary.csv
```

## Yorumlama kuralları

- Repeated spatial CV'de yalnız ortalama değil, MAE %5–%95 aralığı raporlanır.
- Leave-one-domain-out sonuçları BSA içi yapısal extrapolation testidir; insan proteomu
  dış doğrulaması değildir.
- Conditional permutation, sistein komşuluğunun korelasyonlu yapısal özelliklerden
  bağımsız katkısını sınar; nedensellik kanıtlamaz.
- Conformal coverage hedefi nominal olarak %90'dır. Küçük örneklem nedeniyle sonuçlar
  site bazında değil, genel kalibrasyon göstergesi olarak yorumlanır.
- OOD skoru yüksek adaylarda tahmin güveni düşürülmelidir.
