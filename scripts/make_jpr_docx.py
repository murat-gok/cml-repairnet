#!/usr/bin/env python3
"""Convert the achemso manuscript to a Word file for ACS Paragon Plus.

The platform takes DOC, DOCX or ZIP as the manuscript file, not PDF.  This
produces the DOCX from the same LaTeX source that produces the PDF, so the
two cannot drift apart.

achemso defines commands pandoc does not know (tocentry, suppinfo,
acknowledgement, affiliation), and a DOCX cannot embed a PDF image, so the
source is preprocessed before conversion:

* achemso-only environments become ordinary starred sections, in the order
  ACS expects them;
* figures are swapped from vector PDF to 600 dpi PNG;
* the TOC graphic moves to the end under its own heading;
* citations are resolved by citeproc against the ACS CSL style, so the
  reference list matches the compiled PDF.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE.parent / "jpr_tex" / "manuscript_jpr.tex"
BIB = HERE.parent / "jpr_tex" / "references_jpr.bib"
CSL = HERE / "acs.csl"
STAGED = HERE / "manuscript_for_pandoc.tex"
MD = HERE / "manuscript_for_pandoc.md"
OUT = HERE / "manuscript_jpr.docx"

TITLE = (
    "Condition identity and label weighting govern structure-response "
    "inference from reused site-resolved PTM data: a reanalysis of "
    "enzymatic carboxymethyllysine repair"
)


# LaTeX math becomes OMML, which Word renders but other Word-compatible
# readers silently drop -- LibreOffice loses the N of $N^{\varepsilon}$
# outright, leaving "reverses -carboxymethyllysine" on the page.  Every
# expression in this manuscript is simple enough to set as formatted text,
# which renders identically everywhere.  The LaTeX source keeps real math.
MATH_AS_TEXT = {
    r"$N^{\varepsilon}$": r"N\textsuperscript{ε}",
    r"$\mathrm{var}(A)/C^{2} + A^{2}\,\mathrm{var}(C)/C^{4}$":
        "var(A)/C² + A²·var(C)/C⁴",
    r"$(\sum w)^2/\sum w^2$": "(Σw)²/Σw²",
    r"$w_{\min}$": r"w\textsubscript{min}",
    r"$w_{\max}$": r"w\textsubscript{max}",
    r"$\rho = 0.400$": "ρ = 0.400",
    r"$\rho = -0.224$": "ρ = −0.224",
    r"$\alpha = 0.05$": "α = 0.05",
    r"$-0.224$": "−0.224",
    r"$-0.135$": "−0.135",
    r"$\rho$": "ρ",
    r"$\alpha$": "α",
    r"$\Delta$": "Δ",
    r"$>90$": ">90",
    r"$=1$": "=1",
    r"$=0$": "=0",
    r"$P$": r"\textit{P}",
    r"$n$": r"\textit{n}",
}


def _balanced(text: str, start: int) -> tuple[str, int]:
    """Return the content of a brace group starting at ``text[start] == '{'``."""
    depth, i = 0, start
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1 : i], i + 1
        i += 1
    raise ValueError("unbalanced braces")


def number_floats(text: str) -> str:
    """Give tables and figures explicit numbered captions.

    pandoc emits a LaTeX ``\\caption`` both as a paragraph and as the table's
    own caption, and the second copy keeps raw math, so the caption is pulled
    out and rewritten as a single bold paragraph -- above a table, below a
    figure, as ACS sets them.
    """
    counters = {"table": 0, "figure": 0}
    out, cursor = [], 0
    pattern = re.compile(r"\\begin\{(table|figure)\}(\[[^\]]*\])?")
    while True:
        match = pattern.search(text, cursor)
        if match is None:
            out.append(text[cursor:])
            break
        kind = match.group(1)
        end = text.find(r"\end{%s}" % kind, match.end())
        if end == -1:
            out.append(text[cursor:])
            break
        end += len(r"\end{%s}" % kind)
        block = text[match.end() : end]

        caption = ""
        cap_at = block.find(r"\caption")
        if cap_at != -1:
            brace = block.find("{", cap_at)
            caption, after = _balanced(block, brace)
            block = block[:cap_at] + block[after:]
        block = re.sub(r"\\label\{[^}]*\}", "", block)
        block = block.replace(r"\end{%s}" % kind, "")

        counters[kind] += 1
        label = r"\textbf{%s %d.} " % (kind.capitalize(), counters[kind]) + caption

        out.append(text[cursor : match.start()])
        if kind == "table":
            out.append("\n\n" + label + "\n\n" + block + "\n\n")
        else:
            out.append("\n\n" + block + "\n\n" + label + "\n\n")
        cursor = end
    return "".join(out)


def preprocess(source: str) -> str:
    text = source

    # achemso is unknown to pandoc; article carries the same body commands.
    text = re.sub(r"\\documentclass\[[^]]*\]\{achemso\}", r"\\documentclass{article}", text)
    text = text.replace(r"\SectionNumbersOn", "").replace(r"\SectionNumbersOff", "")
    text = text.replace(r"\sloppy", "")

    # The TOC graphic is a separate upload; keep a copy at the end of the
    # document under its own heading, which is where ACS expects it in a
    # Word manuscript, and drop the achemso environment.
    text = re.sub(
        r"\\begin\{tocentry\}.*?\\end\{tocentry\}", "", text, flags=re.DOTALL
    )

    # Front matter.  pandoc's LaTeX reader moves \title, \author and the
    # abstract environment into document metadata, and the markdown stage
    # then drops them -- which silently costs the manuscript its title and
    # its entire abstract.  So the front matter is written into the body.
    affiliation = re.search(r"\\affiliation\{(.+?)\}", text, flags=re.DOTALL)
    email = re.search(r"\\email\{(.+?)\}", text)
    keywords = re.search(r"\\keywords\{(.+?)\}", text, flags=re.DOTALL)
    for command in ("affiliation", "email", "keywords", "abbreviations", "phone", "fax"):
        text = re.sub(r"\\%s\{.*?\}" % command, "", text, flags=re.DOTALL)
    text = re.sub(r"\\title\{.*?\}", "", text, flags=re.DOTALL)
    text = re.sub(r"\\author\{.*?\}", "", text, flags=re.DOTALL)

    # The title and author come from document metadata, which pandoc renders
    # as Word's own Title style; repeating them here would print them twice.
    front = []
    if affiliation:
        front += [affiliation.group(1).strip(), ""]
    if email:
        front += ["Corresponding author: " + email.group(1).strip(), ""]
    front += ["ORCID: 0000-0003-2261-9288", ""]
    text = text.replace(r"\maketitle", "\n".join(front))

    # Keep the abstract as a real section so it survives the markdown stage.
    text = text.replace(r"\begin{abstract}", r"\section*{Abstract}")
    text = text.replace(r"\end{abstract}", "")

    if keywords:
        text = re.sub(
            r"(\\section\*\{Introduction\}|\\section\{Introduction\})",
            r"\\textbf{Keywords:} " + keywords.group(1).strip() + "\n\n" + r"\1",
            text,
            count=1,
        )

    # achemso-only environments -> plain starred sections.
    text = text.replace(r"\begin{suppinfo}", r"\section*{Supporting Information}")
    text = text.replace(r"\end{suppinfo}", "")
    text = text.replace(r"\begin{acknowledgement}", r"\section*{Acknowledgement}")
    text = text.replace(r"\end{acknowledgement}", "")

    # DOCX cannot embed a PDF image.
    text = re.sub(
        r"\{figures/([A-Za-z0-9_]+)\.pdf\}", r"{figures/\1.png}", text
    )

    # citeproc inserts the reference list; the LaTeX command would be copied
    # through as a stray paragraph.
    text = re.sub(r"\\bibliography\{[^}]*\}", "", text)

    # Table 3 groups six of its columns under two spanning headings.  pandoc
    # does not carry \multicolumn across, so the heading lands in one narrow
    # cell and wraps a letter at a time, and the tabular's spare eleventh
    # column becomes a visible empty one.  For Word the grouping is folded
    # into the column names instead, which the caption already explains.
    text = text.replace(r"\begin{tabular}{lrrrrrrrrrr}",
                        r"\begin{tabular}{lrrrrrrrrr}")
    text = re.sub(
        r"[^\n]*\\multicolumn\{3\}\{c\}\{Unweighted \(primary\)\}.*?\\\\\n"
        r"(?:\\cmidrule[^\n]*\n)?"
        r"Load \(ng\)[^\n]*\\\\",
        r"Load (ng) & $n$ & Baseline & MAE (unw.) & $\\rho$ (unw.) & $P$ (unw.)"
        r" & MAE (i.v.) & $\\rho$ (i.v.) & $P$ (i.v.) & SD (IQR) \\\\",
        text,
        flags=re.DOTALL,
        count=1,
    )

    # Set every expression as formatted text rather than OMML.
    for latex_math, replacement in MATH_AS_TEXT.items():
        text = text.replace(latex_math, replacement)

    # pandoc copies \cmidrule's span argument into the first header cell.
    text = re.sub(r"\\cmidrule(\(l?r?\))?\{[0-9]+-[0-9]+\}", "", text)

    # Tables and figures need explicit numbers in a Word manuscript.
    text = number_floats(text)

    # Append the TOC graphic at the end, as ACS asks for in a Word file.
    text = text.replace(
        r"\end{document}",
        "\n".join(
            [
                r"\section*{Table of Contents Graphic}",
                r"\begin{figure}[htbp]\centering",
                r"\includegraphics[width=0.6\textwidth]{figures/toc_graphic.png}",
                r"\end{figure}",
                r"\end{document}",
            ]
        ),
    )
    return text


def main() -> int:
    if not TEX.exists():
        print(f"source not found: {TEX}", file=sys.stderr)
        return 1
    STAGED.write_text(preprocess(TEX.read_text(encoding="utf-8")), encoding="utf-8")

    # Two stages, because citeproc appends the reference list to the very end
    # of the document, which would put it after the TOC graphic.  Going via
    # markdown lets the bibliography be placed explicitly, under its own
    # heading, where ACS expects it.
    markdown = subprocess.run(
        [
            "pandoc", str(STAGED), "--from=latex", "--to=markdown",
            f"--resource-path={HERE}", "--wrap=none",
        ],
        capture_output=True, text=True,
    )
    if markdown.returncode != 0:
        print(markdown.stderr.strip(), file=sys.stderr)
        return markdown.returncode

    body = markdown.stdout
    # pandoc gives every figure the alt text "image", which its docx writer
    # then prints as a line of its own under the picture.
    body = body.replace("![image](", "![](")
    refs_block = "\n\n# References\n\n::: {#refs}\n:::\n\n"
    marker = "# Table of Contents Graphic"
    if marker in body:
        body = body.replace(marker, refs_block.strip() + "\n\n" + marker, 1)
    else:
        body += refs_block
    MD.write_text(body, encoding="utf-8")

    command = [
        "pandoc",
        str(MD),
        "--from=markdown",
        "--to=docx",
        "--standalone",
        "--citeproc",
        f"--bibliography={BIB}",
        f"--csl={CSL}",
        f"--resource-path={HERE}",
        f"--metadata=title:{TITLE}",
        "--metadata=author:Murat Gök",
        f"--output={OUT}",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)
    if result.returncode != 0:
        return result.returncode

    print(f"wrote {OUT} ({OUT.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
