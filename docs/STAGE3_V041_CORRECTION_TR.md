# Stage 3 v0.4.1 Düzeltmesi

## Neden gerekli?

Stage 3 v0.4.0 quick analizinde beş farklı K-means seed'i, küme etiketleri dışında
aynı uzamsal partition'ı üretmiştir. Bunun nedeni `n_init=50` ile her seed'in aynı
global K-means çözümüne yakınsamasıdır. Ayrıca estimator `random_state` değeri
partition seed'iyle birlikte değiştiği için görünen küçük metrik değişimleri yalnız
partition duyarlılığını temsil etmemiştir.

## Düzeltmeler

- K-means tek random initialization (`n_init=1`, `init='random'`) kullanır.
- Partition'lar cluster-label-invariant SHA-256 fingerprint ile karşılaştırılır.
- Duplicate partition'lar atlanır; istenen sayıda benzersiz partition üretilir.
- Model random seed'i bütün partition'larda sabit tutulur.
- Pairwise Adjusted Rand Index (ARI) çeşitlilik raporu yazılır.
- Partition manifesti, grup boyutları ve gerçek kullanılan seed'ler kaydedilir.
- Domain bazlı model metrikleri ve uncertainty özeti ayrıca yazılır.

## Yeni çıktılar

- `results/stage3/repeated_spatial/repeated_spatial_partition_manifest.csv`
- `results/stage3/repeated_spatial/repeated_spatial_partition_diversity.csv`
- `results/stage3/domain_holdout/domain_metrics_by_domain.csv`
- `results/stage3/domain_holdout/domain_uncertainty_by_domain.csv`

## Quick doğrulama

```bash
python -u scripts/run_stage3.py \
  --mode quick \
  --skip-domain \
  --skip-cysteine \
  --n-jobs 1 \
  2>&1 | tee results/stage3_repeated_v041_quick_console.txt
```

Kontrol:

```bash
python - <<'PY'
import pandas as pd
p = "results/stage3/repeated_spatial/repeated_spatial_partition_diversity.csv"
df = pd.read_csv(p)
print(df.to_string(index=False))
assert (df["unique_fingerprints"] == df["partition_count_diagnostic"]).all()
assert (df["pairwise_ari_min"] < 1.0).all()
print("PASS: spatial partitions are genuinely distinct")
PY
```
