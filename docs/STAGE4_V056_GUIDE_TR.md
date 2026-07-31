# Stage 4 v0.5.6 — Pareto ve Domain-Shift Düzeltmesi

Bu yama model tahminlerini veya yapısal özellikleri değiştirmez. Yalnızca önceliklendirme ve audit katmanını düzeltir.

## Değişiklikler

- Pareto sıralaması yalnız `in_domain` adaylarda çalışır.
- Yedi ayrı amaç yerine üç bağımsız amaç kullanılır:
  1. Extra Trees onarılabilirlik tahmini,
  2. deneysel kanıt skoru,
  3. transfer güvenilirliği.
- Transfer güvenilirliği; OOD güvenliği, conformal aralık kesinliği ve model anlaşmasından oluşur.
- Borderline/OOD adaylar ayrı `cml_repairability_domain_shift_watchlist.csv` dosyasına yazılır.
- Audit, domain-shift watchlist sayısını ve in-domain Pareto oranını raporlar.
- Model kartı sürümü `0.5.6` olur.

## Çalıştırma

```bash
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
python scripts/run_stage4.py score --refit-model
python scripts/run_stage4.py prioritize --top-n 50 --max-sites-per-protein 2
python scripts/run_stage4.py audit
```

Yapı ve kaynak indirme aşamalarını yeniden çalıştırmak gerekmez.
