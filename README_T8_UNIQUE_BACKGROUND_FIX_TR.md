# CML-RepairNet T8 benzersiz arka plan düzeltmesi v0.11.2

Bu yama T8 enrichment sağlamlık analizinde aynı gen listesini içeren farklı arka plan
etiketlerinin bağımsız kanıtlar gibi sayılmasını önler.

## Değişiklikler

- Arka plan gen kümeleri SHA-256 ile parmak izlenir.
- Aynı gen kümesine sahip `candidate_universe`, `structure_eligible` veya `matched`
  etiketleri tek benzersiz arka plan olarak sayılır.
- `background_robust`, en az iki benzersiz arka plan kümesinde ve matched analizinde
  anlamlılık gerektirir.
- CORUM/KEGG/Reactome/WikiPathways kök düğümleri `ontology_root_term` olarak işaretlenir.
- `background_robust_informative` yalnız biyolojik olarak yorumlanabilir sağlam terimleri tutar.
- T5–T7, model, atlas skorları ve ağ analizi değiştirilmez.

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Phase2_T8_Unique_Background_Fix_v0.11.2.zip -d .
python -m pytest -q tests/test_revision_phase2.py
```

Beklenen: `6 passed`.

## T8'i yeniden çalıştırma

Mevcut `results/revision_phase2/gprofiler_cache` klasörünü silmeyin.

```bash
set -o pipefail
python -u scripts/run_revision_phase2.py enrichment \
  --controls-per-query 5 \
  2>&1 | tee results/revision_phase2_t8_v0112_console.txt

echo "EXIT_CODE=${PIPESTATUS[0]}"
```

Mevcut cache kullanılacağı için g:Profiler sorgularının yeniden indirilmesi beklenmez.
