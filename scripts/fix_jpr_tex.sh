#!/bin/sh
# md2tex.py ciktisini ACS/kimya gosterimine uygun hale getirir
T="$1"
sed -i 's|\\\\varepsilon|\\varepsilon|g; s|\\\\mathrm|\\mathrm|g; s|\\\\,|\\,|g' "$T"
sed -i 's|C-alpha|C$\\alpha$|g' "$T"
sed -i 's|^\\usepackage{microtype}|\\usepackage{lmodern}|' "$T"
grep -q '^\\sloppy' "$T" || sed -i 's|^\\begin{document}|\\sloppy\n\\begin{document}|' "$T"
