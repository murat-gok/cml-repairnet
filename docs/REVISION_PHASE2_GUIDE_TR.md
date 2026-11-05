# CML-RepairNet Revision Phase II — T5–T8

Bu paket, PEDS hedefli nihai bilimsel dondurmadan önce dört duyarlılık analizi ekler.

## Bilimsel sınır

Paket yeni site-resolved deney etiketi üretmez ve modeli dış sonuçlara göre yeniden seçmez.
Tüm atlas aralıkları ve işlevsel analizler hipotez üretme amaçlıdır.

## T5 — Covariate-shift ve conformal duyarlılığı

Karşılaştırılan aralıklar:

- standart OOF residual aralığı,
- atlas/BSA dağılım ayrımından hesaplanan yaklaşık density-ratio ağırlıklı aralık,
- yerel yapısal uzaklık ağırlıklı aralık,
- yalnız betimleyici OOD-inflated aralık.

Density-ratio ve yerel aralıklar duyarlılık analizidir; bağımsız proteinlerde finite-sample coverage garantisi olarak sunulmaz.

## T6 — Boundary shrinkage

Altı uzamsal partisyonun site-bazlı ortalama OOF tahminleri kullanılır. Kalibrasyon eğimi,
düşük/yüksek yanıt bölgelerinde yönlü hata, leave-one-site-out affine/isotonic diagnostikleri ve
site bootstrap güven aralıkları yazılır. Recalibrated değerler atlası yeniden eğitmek için kullanılmaz.

## T7 — Kompartıman ve erişilebilirlik

UniProt subcellular-location metinleri çoklu bayraklarla ayrıştırılır. Extracellular/surface,
secretory-lumen, intracellular, mitochondrial ve mixed/unknown katmanları oluşturulur. Hiçbir aday
salt kompartıman nedeniyle silinmez; intracellular adaylar secondary hypothesis olarak işaretlenir.

## T8 — Çoklu arka plan enrichment ve STRING çalışma yanlılığı

g:Profiler her sorguyu dört özel arka planla çalıştırır:

1. bütün yayımlanmış aday evreni,
2. yapı-uygun proteinler,
3. in-domain proteinler,
4. protein uzunluğu, lizin sayısı, kanıt ve yapı kalitesiyle eşleştirilmiş arka plan.

Ayrıca STRING ağları text-mining-dominated kenarlar çıkarılarak yeniden özetlenir.

## Kurulum

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Revision_Phase2_v0.11.0.zip -d .
python -m pip install -e ".[dev,experiment]"
python -m pytest -q tests/test_revision_phase2.py
# İsteğe bağlı tam regresyon: python -m pytest -q
```

## Ön kontrol

```bash
python scripts/run_revision_phase2.py doctor
```

## Önerilen sıra

Önce yalnız T5 çalıştırılır; çıktı incelendikten sonra T6–T8’e geçilir.

```bash
mkdir -p results/revision_phase2
set -o pipefail

python -u scripts/run_revision_phase2.py conformal \
  --alpha 0.10 \
  2>&1 | tee results/revision_phase2_t5_console.txt

python -u scripts/run_revision_phase2.py boundary \
  --bootstrap-iterations 20000 \
  2>&1 | tee results/revision_phase2_t6_console.txt

python -u scripts/run_revision_phase2.py compartments \
  2>&1 | tee results/revision_phase2_t7_console.txt

python -u scripts/run_revision_phase2.py enrichment \
  --controls-per-query 5 \
  2>&1 | tee results/revision_phase2_t8_console.txt
```

T8 internet bağlantısı gerektirir ve g:Profiler yanıtlarını `results/revision_phase2/gprofiler_cache/`
altında saklar. Bağlantı olmadan yalnız manifest ve ağ duyarlılığı üretmek için `--offline` kullanılır.
