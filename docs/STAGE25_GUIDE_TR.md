# CML-RepairNet Stage 2.5 Kullanım Rehberi

## Amaç

Stage 2.5, Stage 2'de saptanan yapı-ağırlıklı sinyali aşağıdaki ek denetimlerle sınar:

- Güçlü düzenlileştirilmiş XGBoost ve LightGBM
- Extra Trees ve Ridge ana referansları; Random Forest ve Elastic Net isteğe bağlı kontroller
- Outer-fold dışı etiket kullanmayan Ridge–Extra Trees convex ensemble
- Ölçüm belirsizliğine dayalı örnek ağırlığı duyarlılığı
- Eğitim fold'u içinde korelasyon tabanlı özellik azaltımı
- Held-out permutation importance
- İsteğe bağlı SHAP açıklamaları

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

## Ön kayıtlı birincil analiz

```bash
python -u scripts/run_stage25.py --mode quick 2>&1 \
  | tee results/stage25_quick_console.txt
```

Birincil koşul sabittir:

- Hedef: `repair_published`
- Uzamsal grup: 4
- Örnek ağırlığı: yok
- Özellik paneli: 11 yapısal özellik

Bu analiz, XGBoost veya LightGBM görüldükten sonra koşul değiştirilmesini önler.

## Tam sağlamlık analizi

```bash
set -o pipefail
python -u scripts/run_stage25.py \
  --mode full \
  --bootstrap 20000 \
  --permutations 50000 \
  --permutation-repeats 100 \
  --n-jobs 1 \
  2>&1 | tee results/stage25_full_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

Tam mod şunları çalıştırır:

1. İki hedef × 3/4/5 uzamsal grup, tam yapı paneli, ağırlıksız.
2. Birincil koşulda clipped ve sqrt örnek ağırlıkları.
3. Birincil koşulda nested korelasyon azaltımı.

Tüm kombinasyonların çapraz ürünü özellikle çalıştırılmaz. Bu, 33 örnekte model-selection bias ve gereksiz hesaplama yükünü sınırlar.

## SHAP

Önce Stage 2.5 sonucunda en iyi tree modelini belirleyin. Ardından örneğin:

```bash
python scripts/run_stage25_shap.py --model extra_trees
python scripts/run_stage25_shap.py --model xgboost
python scripts/run_stage25_shap.py --model lightgbm
```

SHAP çıktısı performans doğrulaması değildir; tüm veri üzerinde kurulan nihai modelin betimleyici açıklamasıdır.

## Başlıca çıktılar

```text
results/stage25/stage25_all_metrics.csv
results/stage25/stage25_all_oof_predictions.csv
results/stage25/stage25_paired_comparisons.csv
results/stage25/stage25_reference_comparisons.csv
results/stage25/stage25_primary_permutation_importance.csv
results/stage25/stage25_report.md
```

Her koşulun ayrıntıları kendi alt klasöründe saklanır. `stage25_fold_details.json`, her dış fold için seçilen hiperparametreleri, özellikleri ve ensemble ağırlığını içerir.

## Bilimsel karar kuralları

Bir boosting modelini ana modele yükseltmek için:

- Birincil koşulda Extra Trees/Ridge'den daha düşük OOF MAE,
- Fark için destekleyici paired bootstrap ve sign-flip sonucu,
- İki hedef ve 3/4/5 grup düzenlerinde aynı yön,
- Clipped/sqrt ağırlıklarda performans çökmesi olmaması,
- Önem sıralamasının biyolojik açıdan savunulabilir olması

gerekir.

Tek bir en düşük skor yeterli değildir.

## İsteğe bağlı geniş model seti

Varsayılan koşu küçük veri için ana adaylara odaklanır: dummy, Ridge, Extra Trees, XGBoost ve LightGBM. Random Forest ile Elastic Net'i de yeniden çalıştırmak için:

```bash
python -u scripts/run_stage25.py --mode quick --models \
  dummy_mean ridge elastic_net random_forest extra_trees xgboost lightgbm
```
