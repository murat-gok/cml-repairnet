# Stage 6 v0.7.1 düzeltme rehberi

## Neden gerekli?

Stage 4 portföyü tür dengesi ve protein başına en fazla iki bölge kısıtıyla
oluşturuldu. Stage 6 v0.7.0 ise stres iterasyonlarında kısıtsız global ilk 100'ü
seçti ve Pareto cephesini insan/fare birlikte hesapladı. Bu nedenle
`top100_selection_frequency`, gerçek Stage 4 portföy kararlılığı değildi ve Tier
etiketlerini etkiliyordu.

v0.7.1:

- Pareto cephesini her tür içinde ayrı hesaplar.
- Her iterasyonda Stage 4 tür kotasını ve protein başına iki bölge sınırını
  yeniden uygular.
- `portfolio_selection_frequency` üretir.
- Tier 1/2 sınıflarını bu yeni frekansa göre belirler.
- Kısıtsız top-k frekanslarını yalnız tanısal çıktı olarak korur.

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate

unzip -o \
  /mnt/c/Users/user/Downloads/CML_RepairNet_Stage6_Portfolio_Stability_Fix_v0.7.1.zip \
  -d .

python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

## Yeniden çalıştırılacak aşamalar

Conservation, enrichment ve network sonuçlarını tekrar üretmek gerekmez.

```bash
python -u scripts/run_stage6.py rank \
  2>&1 | tee results/stage6_rank_v071_console.txt

python -u scripts/run_stage6.py tiers \
  2>&1 | tee results/stage6_tiers_v071_console.txt

python scripts/run_stage6.py figures
python scripts/run_stage6.py report
```

Alternatif olarak tüm Stage 6'yı yeniden çalıştırabilirsiniz:

```bash
set -o pipefail
python -u scripts/run_stage6.py all \
  2>&1 | tee results/stage6_v071_full_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

## Hızlı doğrulama

```bash
python - <<'PY'
import json
import pandas as pd

rank = pd.read_csv(
    "data/processed/stage6/stage6_rank_stability.csv"
)
tiers = pd.read_csv(
    "data/processed/stage6/stage6_candidate_tiers.csv"
)
with open("results/stage6/stage6_audit.json") as f:
    audit = json.load(f)

assert "portfolio_selection_frequency" in rank.columns
assert audit["stage6_version"] == "0.7.1"

expected_slots = rank.groupby("species")[
    "portfolio_selection_frequency"
].sum()
print("İterasyon başına beklenen portföy bölgesi:")
print(expected_slots.to_string())

print("\nTier sayıları:")
print(tiers["candidate_tier"].value_counts().to_string())

print("\nPASS: Stage 6 portfolio stability is Stage-4-compatible")
PY
```

Aynı dondurulmuş girdilerle beklenen frekans toplamları:

```text
human    50.0
mouse    50.0
```

Referans doğrulama çalışmasında Tier dağılımı:

```text
Tier 1: 4
Tier 2: 33
Tier 3: 712
Watchlist: 519
```
