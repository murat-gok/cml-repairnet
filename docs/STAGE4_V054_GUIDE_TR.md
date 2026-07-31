# Stage 4 v0.5.4 Kullanım Rehberi

## Uygulama

```bash
cd /mnt/c/Users/user/RepairNet
source .venv/bin/activate
unzip -o /mnt/c/Users/user/Downloads/CML_RepairNet_Stage4_Parser_Evidence_Fix_v0.5.4.zip -d .
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
```

Kullanıcının v0.5.3 kurulumunda 33 test geçtiği için bu yama sonrasında 39 test
beklenir.

## Eski Stage 4 ara çıktılarını temizleme

```bash
rm -f data/interim/stage4/stage4_candidates.csv
rm -f data/interim/stage4/stage4_prm_evidence.csv
rm -f data/interim/stage4/stage4_candidates_uniprot.csv
rm -f data/interim/stage4/stage4_prm_mapping.csv
rm -f results/stage4/stage4_parse_diagnostics.csv
rm -f results/stage4/stage4_excluded_sequence.csv
rm -f results/stage4/stage4_excluded_prm_mapping.csv
```

UniProt ve AlphaFold cache klasörleri silinmemelidir.

## Quick prepare

```bash
python -u scripts/run_stage4.py prepare --mode quick \
  2>&1 | tee results/stage4_prepare_v054_console.txt
```

Ham dosyaların doğrulanmış tam parse kapsamı:

- MEF discovery sites: 1092
- HUVEC discovery sites: 298
- Human/mouse overlap evidence rows: 302
- Heart aging sites: 180
- Kidney aging sites: 59
- Liver aging sites: 94
- Tissue discovery sites: 374
- Proteostasis sites: 7
- Mouse PRM peptides with CML marker: 16
- HUVEC PRM peptides with CML marker: 1

Tam tabloda 1681 benzersiz aday ve 778 protein bulunmuştur. Quick seçim 25 protein
içinde 54 site üretir; PRM ile haritalanan yeni siteler annotate aşamasında ayrıca
eklenebilir.

## Annotate

```bash
python -u scripts/run_stage4.py annotate \
  2>&1 | tee results/stage4_annotate_v054_console.txt
```

Yeni denetim dosyaları:

- `data/interim/stage4/stage4_prm_evidence.csv`
- `data/interim/stage4/stage4_prm_mapping.csv`
- `results/stage4/stage4_excluded_prm_mapping.csv`

PRM peptitleri canonical UniProt dizisinde aranır; CML işaretli lizinin mutlak
konumu buradan hesaplanır. Exact eşleme bulunmazsa I/L eşdeğerliği denenir.
Çoklu eşleme veya dizi uyuşmazlığı açık gerekçeyle dışlanır.
