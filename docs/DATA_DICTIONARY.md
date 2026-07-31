# Data dictionary

## Primary unit of analysis

One experimentally observed **CML-modified lysine site**.

## Site numbering

The article's Figure 3 labels BSA sites in the mature-protein numbering system. UniProt
P02769 contains a 24-residue signal/propeptide segment before the mature chain. Therefore:

```text
site_uniprot = site_paper + 24
```

Both identifiers are retained. Every mapped UniProt residue must validate as `K`.

## Primary label

- `repair_published`: site-level repair fraction from article Figure 3C, range 0–1.

## Replicate-derived labels

- `active_rep1..3`: percentage CML abundance after CMLase treatment.
- `control_rep1..3`: percentage CML abundance after inactive-enzyme treatment.
- `repair_from_means_raw = 1 - mean(active)/mean(control)`.
- `repair_from_means_clipped`: raw value clipped to 0–1.
- `repair_variance_approx`: delta-method approximation using replicate variability.
- `label_weight`: inverse-variance weight normalized to median 1.

The published value remains the primary endpoint. Raw and replicate-derived values are used for
sensitivity analyses and heteroscedastic modelling.

## Core feature groups

1. Local sequence composition and physicochemical context.
2. Frozen protein-language-model residue embedding.
3. Experimental structure: B-factor, solvent-accessible area, packing/contact counts.
4. Spatial coordinates used only to form blocked validation groups.
5. Later phases: pLDDT/PAE, electrostatics, residue depth, ECM/tissue annotations and docking-derived accessibility.
