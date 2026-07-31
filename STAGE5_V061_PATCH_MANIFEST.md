# CML-RepairNet Stage 5 v0.6.1 Patch

- Corrects g:Orth field mapping (`name` target gene symbol, `converted` source ENSG,
  `ortholog_ensg` target ENSG).
- Adds case-insensitive human/mouse symbol matching and duplicate-pair removal.
- Prevents zero-byte cross-species summary files.
- Adds orthology schema diagnostics to the manifest and audit.
- Adds module silhouette quality to the audit and marks low-separation modules as exploratory.
- Adds species/query/source context to enrichment figure labels.
- Adds two regression tests for the current official g:Orth schema.
