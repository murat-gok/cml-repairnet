# CML-RepairNet — VS Code Kurulum ve İlk Çalıştırma

## 1. Önerilen çalışma ortamı

Windows üzerinde **VS Code + Remote WSL + Python 3.11** kullanın. Projeyi WSL'nin Linux dosya
sisteminde tutun:

```bash
mkdir -p ~/projects
cd ~/projects
unzip /mnt/c/Users/user/Downloads/CML-RepairNet_VSCode.zip
cd CML-RepairNet
code .
```

`/mnt/c/...` altında çalışmak mümkündür; ancak çok sayıda küçük dosya, model önbelleği ve ilerideki
MD/PLM iş yüklerinde Linux dosya sistemi genellikle daha uygundur.

## 2. VS Code eklentileri

Depo açıldığında önerilen eklentileri yükleyin:

- Python
- Pylance
- Ruff
- Jupyter
- Remote - WSL

Bunlar `.vscode/extensions.json` içinde tanımlıdır.

## 3. Ortamı oluşturma

VS Code menüsünden:

**Terminal → Run Task → CML: 0 - Bootstrap environment**

veya terminalden:

```bash
bash scripts/bootstrap.sh
```

Komut şu işlemleri yapar:

1. `.venv` Python 3.11 ortamını oluşturur.
2. Projeyi editable modda kurar.
3. Geliştirme/test bağımlılıklarını yükler.
4. Birim testlerini çalıştırır.

Python 3.11 bulunamazsa Ubuntu/WSL dağıtımınıza uygun şekilde Python 3.11 ve `venv` desteğini
kurmanız gerekir. Ardından aynı komutu yeniden çalıştırın.

## 4. Python yorumlayıcısını seçme

VS Code'da `Ctrl+Shift+P`:

```text
Python: Select Interpreter
```

Şunu seçin:

```text
<proje>/.venv/bin/python
```

WSL dışında doğrudan Windows ile çalışıyorsanız yol `.venv\Scripts\python.exe` olur; ancak proje
için WSL akışı önerilir.

## 5. Ortam denetimi

```bash
source .venv/bin/activate
cml-repairnet doctor
```

İlk çalıştırmada ham veriler bulunmadığı için `NOT READY FOR STAGE 1` görülmesi normaldir.

## 6. İlk bilimsel aşama: veri temelini kurma

```bash
cml-repairnet download
cml-repairnet stage1 --no-download
```

İlk komut yalnızca resmi açık kaynaklardan şu girdileri indirir:

- makale kaynak XLSX dosyası,
- UniProt P02769 FASTA dizisi,
- RCSB 4F5S PDB yapısı.

İkinci komut:

- 33 BSA bölgesini çıkarır,
- UniProt konumlarının lizin olduğunu doğrular,
- makaledeki olgun BSA numaralarını 4F5S yapısına eşler,
- dizi ve yapı özelliklerini üretir,
- veri bütünlüğü raporunu yazar.

Başarılı sonuç:

```text
Stage 1 completed and audited.
Sites: 33
Report: results/stage1_data_audit.md
```

## 7. VS Code görevleri

`Terminal → Run Task` altında:

- `CML: Doctor`
- `CML: 1 - Download public data`
- `CML: 2 - Run Stage 1 data foundation`
- `CML: Tests`
- `CML: Train leakage-controlled baselines`

bulunur.

## 8. Hata ayıklama

Run and Debug panelinde:

- `CML: Doctor`
- `CML: Stage 1 (no download)`
- `Python: Current file`

profilleri hazırdır. Bir satırın soluna tıklayarak breakpoint koyabilirsiniz.

## 9. İlk kilometre taşı tamamlanma ölçütleri

Stage 1 tamamlanmış sayılmak için:

- tam 33 etiketli bölge bulunmalı,
- 33 UniProt konumunun tamamı `K` olmalı,
- 4F5S eşlemesinde eksik bölge olmamalı,
- hedefler `[0,1]` aralığında olmalı,
- üç tablonun site anahtarları bire bir eşleşmeli,
- `results/stage1_data_audit.md` durumu `PASS` olmalı.

Bu koşullar geçmeden ESM-2, XGBoost, docking veya MD aşamasına geçilmez.
