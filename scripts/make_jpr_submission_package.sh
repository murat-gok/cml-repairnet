#!/bin/sh
# Build the three files ACS Paragon Plus asks for.
#
#   Manuscript File*  DOC / DOCX / ZIP  -> manuscript_jpr.docx (preferred)
#                                       or manuscript_jpr.zip (LaTeX source)
#   Cover Letter*     DOCX/DOC/PDF/RTF/TXT -> cover_letter_jpr.pdf
#   Other Files       (Supporting Information for publication) -> supporting_information.pdf
#
# The platform rejects a PDF as the Manuscript File, so the achemso source goes
# up as a ZIP: .tex, .bib, the generated .bbl (so the bibliography renders even
# if BibTeX is not run), and the figures as both vector PDF and 600 dpi TIFF.
set -e

SP="$(cd "$(dirname "$0")" && pwd)"
TEX="$SP/jpr_tex"
OUT="$SP/submission_package"

rm -rf "$OUT" "$SP/_zipstage"
mkdir -p "$OUT" "$SP/_zipstage/figures"

# --- Manuscript File: LaTeX source archive -------------------------------
cp "$TEX/manuscript_jpr.tex"  "$SP/_zipstage/"
cp "$TEX/references_jpr.bib"  "$SP/_zipstage/"
cp "$TEX/manuscript_jpr.bbl"  "$SP/_zipstage/"
for f in toc_graphic figure1_reuse_workflow figure2_label_weighting figure3_replication_composition; do
  cp "$TEX/figures/$f.pdf" "$SP/_zipstage/figures/"
  [ -f "$TEX/figures/$f.tif" ] && cp "$TEX/figures/$f.tif" "$SP/_zipstage/figures/"
done

( cd "$SP/_zipstage" && zip -q -r -X "$OUT/manuscript_jpr.zip" . )
rm -rf "$SP/_zipstage"

# --- Manuscript File (preferred): Word -------------------------------
# Built from the same LaTeX source by docx/make_docx.py, so the Word file
# and the PDF cannot drift apart.
cp "$SP/docx/manuscript_jpr.docx" "$OUT/"

# --- Cover Letter ---------------------------------------------------------
cp "$SP/cl/cover_letter_jpr.pdf" "$OUT/"

# --- Supporting Information ----------------------------------------------
cp "$SP/si/supporting_information.pdf" "$OUT/"

# --- Reference copy of the compiled manuscript (not uploaded) ------------
cp "$TEX/manuscript_jpr.pdf" "$OUT/manuscript_jpr_REFERENCE_ONLY.pdf"

echo "--- $OUT ---"
ls -la "$OUT" | tail -n +4 | awk '{printf "  %9s  %s\n",$5,$9}'
echo
echo "--- manuscript_jpr.zip contents ---"
unzip -l "$OUT/manuscript_jpr.zip" | sed 's/^/  /'
