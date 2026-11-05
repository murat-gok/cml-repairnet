# CML-RepairNet Final Scientific Freeze v0.12.0

Bu aşama yeni model eğitmez ve adayları yeniden seçmez. Frozen Stage 4–6 skor, sıralama, Pareto, portföy ve tier sonuçlarını T5–T8 belirsizlik, kalibrasyon, kompartıman ve fonksiyonel sağlamlık sonuçlarıyla birleştirir.

## Bilimsel sınırlar

- Ana çıktı `cmlase_response_score`: göreli önceliklendirme skorudur; mutlak onarım oranı değildir.
- `bsa_calibrated_sensitivity_*`: BSA içinde kalibre edilmiş duyarlılık bandıdır; proteom transferinde coverage garantisi değildir.
- `ood_inflated_descriptive_*`: OOD arttıkça genişleyen betimsel banttır; formal conformal aralık değildir.
- `candidate_tier`: hesaplamalı öncelik düzeyidir.
- `biomedical_interpretation_class`: kompartıman, erişilebilirlik ve teslimat bağlamını ayrı gösterir.
- External scaffold analizi dış doğrulama değil applicability stress testidir.

## Çalıştırma

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate

python -m pytest -q tests/test_final_scientific_freeze.py
python scripts/run_final_scientific_freeze.py doctor
python scripts/run_final_scientific_freeze.py build
```

## Çıktılar

`results/final_scientific_freeze/` altında:

- `cml_repairnet_final_atlas_all_sites.csv`
- `cml_repairnet_final_atlas_curated.csv`
- `cml_repairnet_final_manuscript_candidates.csv`
- `cml_repairnet_final_priority_portfolio.csv`
- `cml_repairnet_final_pareto_candidates.csv`
- `cml_repairnet_final_domain_shift_watchlist.csv`
- `cml_repairnet_final_background_robust_enrichment.csv`
- `cml_repairnet_final_string_sensitivity.csv`
- `final_freeze_manifest.json`
- `FINAL_SCIENTIFIC_FREEZE_REPORT.md`
