# CML-RepairNet Stage 6 v0.7.2 düzeltme rehberi

## Neden gerekli?

Stage 6 v0.7.1 tür dengesini ve protein başına iki bölge sınırını yeniden
uyguluyordu; ancak yalnızca ilk Pareto cephesini ayırıyor, ikinci ve sonraki
nondominated katmanları tek bir "dominated" havuzunda birleştiriyordu.

Dondurulmuş Stage 4 seçimi hiçbir pertürbasyon olmadan yeniden çalıştırıldığında
bu yaklaşım 100 portföy bölgesinin yalnızca 70'ini geri üretiyordu. Stage 4'ün
gerçek seçim sırası şöyledir:

1. İnsan ve fare için nondominated front rank ayrı hesaplanır.
2. Daha küçük front rank önce gelir (`front 1`, `front 2`, ...).
3. Aynı front içindeki adaylar weighted priority azalan biçimde sıralanır.
4. Tür başına 50 bölge ve protein başına en fazla iki bölge sınırı uygulanır.

v0.7.2 bu tam sıralamayı her stres iterasyonunda yeniden hesaplar ve frozen
Stage 4 portföyünü pertürbasyonsuz durumda 100/100 geri üretir.

## Uygulama

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate

unzip -o \
  /mnt/c/Users/user/Downloads/CML_RepairNet_Stage6_FullFront_Stability_Fix_v0.7.2.zip \
  -d .

python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

v0.7.1'de 66 test geçtiyse üç yeni testle yaklaşık 69 test beklenir. Sayının
kendisi değil, bütün testlerin geçmesi esastır.

## Yalnız gerekli aşamaları yeniden çalıştırma

Conservation, enrichment ve network sonuçları değişmez.

```bash
python -u scripts/run_stage6.py rank \
  2>&1 | tee results/stage6_rank_v072_console.txt

python -u scripts/run_stage6.py tiers \
  2>&1 | tee results/stage6_tiers_v072_console.txt

python scripts/run_stage6.py figures
python scripts/run_stage6.py report
```

## Doğrulama

```bash
python - <<'PY'
import json
import pandas as pd

rank = pd.read_csv("data/processed/stage6/stage6_rank_stability.csv")
tiers = pd.read_csv("data/processed/stage6/stage6_candidate_tiers.csv")
weights = pd.read_csv("data/processed/stage6/stage6_priority_weight_sensitivity.csv")

with open("results/stage6/stage6_audit.json") as handle:
    audit = json.load(handle)

print("Stage 6 version:", audit["stage6_version"])
print("Completion:", audit["stage6_completion_status"])
print("Baseline reproduction:",
      audit["baseline_portfolio_reproduction_count"], "/ 100")
print("\nSpecies capacity:")
print(rank.groupby("species")["portfolio_selection_frequency"].sum())
print("\nTier counts:")
print(tiers["candidate_tier"].value_counts())
print("\nWeight scenarios:")
print(weights[["scenario", "portfolio_overlap_count", "portfolio_jaccard"]])

assert audit["stage6_version"] == "0.7.2"
assert audit["stage6_completion_status"] == "complete"
assert audit["baseline_portfolio_reproduction_count"] == 100
assert audit["baseline_portfolio_reproduction_fraction"] == 1.0
assert rank.groupby("species")["portfolio_selection_frequency"].sum().round(8).to_dict() == {
    "human": 50.0,
    "mouse": 50.0,
}
print("\nPASS: Stage 6 full-front portfolio stability is valid")
PY
```

Aynı frozen girdiler ve seed ile beklenen Tier dağılımı:

```text
Tier 3       698
Watchlist    519
Tier 2        42
Tier 1         9
```

## Nihai sonuç paketi

```bash
zip -r CML_RepairNet_Stage6_v072_final_results.zip \
  data/processed/stage6 \
  results/stage6 \
  results/stage6_rank_v072_console.txt \
  results/stage6_tiers_v072_console.txt \
  configs/stage6.yaml \
  tests/test_stage6.py \
  src/cml_repairnet/robustness \
  scripts/run_stage6.py
```
