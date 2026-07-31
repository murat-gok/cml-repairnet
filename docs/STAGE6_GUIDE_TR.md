# CML-RepairNet Stage 6 v0.7.0 kullanım rehberi

## Amaç

Stage 6, Stage 4 ve Stage 5'in dondurulmuş sonuçlarını sağlamlık açısından
sınar. Yeni model aramaz ve 33 bölgeli BSA eğitim setini yeniden optimize etmez.

## Kurulum

ZIP'i proje köküne açın:

```bash
cd /mnt/c/Users/user/RepairNet
unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Stage6_Patch_v0.7.0.zip -d .
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

Stage 5 v0.6.1 dondurmasında 57 test geçtiyse, eklenen 7 Stage 6 testiyle
birlikte **64 passed** beklenir.

## Girdi kontrolü

```bash
python scripts/run_stage6.py doctor
```

`ready: true` beklenir.

## Tek komutla çalıştırma

```bash
set -o pipefail
python -u scripts/run_stage6.py all 2>&1 | tee results/stage6_full_console.txt
echo "EXIT_CODE=${PIPESTATUS[0]}"
```

Varsayılan 2.000 rank-stress ve 500 Pareto-stress iterasyonu hesaplama
açısından ağırdır. Teknik kurulum kontrolü için `configs/stage6.yaml` içindeki
nihai iterasyon değerlerini değiştirmeyin; bunun yerine aşağıdaki komutları
kullanın:

```bash
python scripts/run_stage6.py doctor
python -m pytest -q tests/test_stage6.py
```

Nihai sonuç üretiminde config değerlerini değiştirmeyin.

## Kontrollü çalıştırma

```bash
python scripts/run_stage6.py rank
python scripts/run_stage6.py conservation
python scripts/run_stage6.py enrichment
python scripts/run_stage6.py network
python scripts/run_stage6.py tiers
python scripts/run_stage6.py figures
python scripts/run_stage6.py report
```

## Ana çıktılar

```text
data/processed/stage6/stage6_rank_stability.csv
data/processed/stage6/stage6_priority_weight_sensitivity.csv
data/processed/stage6/stage6_residue_conservation.csv
data/processed/stage6/stage6_enrichment_driver_sensitivity.csv
data/processed/stage6/stage6_network_threshold_sensitivity.csv
data/processed/stage6/stage6_network_hub_removal.csv
data/processed/stage6/stage6_network_degree_preserving_null.csv
data/processed/stage6/stage6_candidate_tiers.csv
data/processed/stage6/stage6_manuscript_candidate_table.csv
results/stage6/stage6_audit.json
results/stage6/stage6_robustness_report.md
results/stage6/figures/
```

## Değerlendirme ZIP'i

```bash
zip -r CML_RepairNet_Stage6_v070_results.zip \
  data/processed/stage6 \
  results/stage6 \
  results/stage6_full_console.txt \
  configs/stage6.yaml \
  tests/test_stage6.py \
  src/cml_repairnet/robustness \
  scripts/run_stage6.py
```
