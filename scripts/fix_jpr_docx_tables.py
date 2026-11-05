#!/usr/bin/env python3
"""Size the manuscript's table columns to their contents.

pandoc gives every column in a table the same width (1584 twips).  That is
fine for a narrow table, but Table 3 carries ten columns, so the row comes
to about eleven inches, overflows the text block, and Word breaks the
numbers across lines -- "0.337" is set as "0.33" then "7", which is
unreadable and looks like a corrupted file.

This rewrites each table's grid so the columns share the available text
width in proportion to the longest cell in each, with a floor so a narrow
column is still legible.  Widths are also written onto the cells and the
table is put in fixed layout, because Word and LibreOffice disagree about
what to do when only the grid is specified.
"""

from __future__ import annotations

import shutil
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
NS = {"w": W}
ET.register_namespace("w", W)

#: A4 with the default one-inch margins, in twentieths of a point.
TEXT_WIDTH = 9360
#: No column narrower than this, so short headers stay on one line.
MIN_WIDTH = 620
#: Rows treated as headers, whose text may wrap.
HEADER_ROWS = 2
#: Longest header length allowed to influence a column's width.
HEADER_CAP = 7


def tag(name: str) -> str:
    return f"{{{W}}}{name}"


def cell_text(cell: ET.Element) -> str:
    return "".join(node.text or "" for node in cell.iter(tag("t")))


def span(cell: ET.Element) -> int:
    merge = cell.find("w:tcPr/w:gridSpan", NS)
    if merge is None:
        return 1
    return int(merge.get(tag("val"), "1"))


def column_demands(table: ET.Element, n_columns: int) -> list[float]:
    """Longest cell content per column, ignoring spanning cells."""
    demands = [0.0] * n_columns
    for position, row in enumerate(table.findall("w:tr", NS)):
        index = 0
        for cell in row.findall("w:tc", NS):
            width = span(cell)
            if width == 1 and index < n_columns:
                length = len(cell_text(cell).strip())
                # A header may wrap onto a second line without harm, but a
                # number must not, so a long header does not get to claim
                # width that the data in its column does not need.
                if position < HEADER_ROWS:
                    length = min(length, HEADER_CAP)
                demands[index] = max(demands[index], length)
            index += width
    # An empty column still needs room for its rule.
    return [max(d, 1.0) for d in demands]


def resize(table: ET.Element) -> bool:
    grid = table.find("w:tblGrid", NS)
    if grid is None:
        return False
    columns = grid.findall("w:gridCol", NS)
    if len(columns) < 2:
        return False

    demands = column_demands(table, len(columns))
    # Allocate the floor first, then share the rest by demand.
    spare = TEXT_WIDTH - MIN_WIDTH * len(columns)
    if spare < 0:  # more columns than will fit at the floor; split evenly
        widths = [TEXT_WIDTH // len(columns)] * len(columns)
    else:
        total = sum(demands)
        widths = [int(MIN_WIDTH + spare * d / total) for d in demands]
    # Put any rounding remainder in the widest column.
    widths[widths.index(max(widths))] += TEXT_WIDTH - sum(widths)

    for column, width in zip(columns, widths):
        column.set(tag("w"), str(width))

    properties = table.find("w:tblPr", NS)
    if properties is not None:
        for existing in properties.findall("w:tblLayout", NS):
            properties.remove(existing)
        layout = ET.SubElement(properties, tag("tblLayout"))
        layout.set(tag("type"), "fixed")

    # Word trusts the cell widths over the grid, so set those too.
    for row in table.findall("w:tr", NS):
        index = 0
        for cell in row.findall("w:tc", NS):
            width = span(cell)
            total = sum(widths[index : index + width]) or MIN_WIDTH
            properties = cell.find("w:tcPr", NS)
            if properties is None:
                properties = ET.Element(tag("tcPr"))
                cell.insert(0, properties)
            for existing in properties.findall("w:tcW", NS):
                properties.remove(existing)
            cell_width = ET.SubElement(properties, tag("tcW"))
            cell_width.set(tag("w"), str(total))
            cell_width.set(tag("type"), "dxa")
            index += width
    return True


def main(path: str) -> int:
    source = Path(path)
    backup = source.with_suffix(".prefit.docx")
    shutil.copy2(source, backup)

    with zipfile.ZipFile(backup) as archive:
        names = archive.namelist()
        payload = {name: archive.read(name) for name in names}

    root = ET.fromstring(payload["word/document.xml"].decode("utf-8"))
    resized = sum(resize(table) for table in root.iter(tag("tbl")))
    payload["word/document.xml"] = ET.tostring(root, encoding="utf-8", xml_declaration=True)

    with zipfile.ZipFile(source, "w", zipfile.ZIP_DEFLATED) as out:
        for name in names:
            out.writestr(name, payload[name])

    print(f"{resized} table(s) resized to fit the text width")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "manuscript_jpr.docx"))
