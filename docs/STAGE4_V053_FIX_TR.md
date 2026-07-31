# Stage 4 v0.5.3 düzeltmesi

Bu düzeltme, protein adlarındaki 10 harfli büyük kelimelerin yanlışlıkla UniProt
accession kabul edilmesini önler. Prepare çıktısının yeniden oluşturulması gerekir.

Kontrol için `uniprot_accession_raw` sütunundaki tüm değerler resmi UniProtKB 6 veya
10 karakter biçimine uymalıdır. Isoform eki (`-2` gibi) kabul edilir ve API sorgusunda
kanonik accession kullanılır.
