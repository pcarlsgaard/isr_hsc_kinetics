"""Minimal standard-library .xlsx reader (no openpyxl in the lock): returns {sheet_name: [[cell, ...], ...]}."""
import re, zipfile, xml.etree.ElementTree as ET
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
def _col(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0): n = n * 26 + ord(ch) - 64
    return n - 1
def read_xlsx(path):
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            ss.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    out = {}
    for sh in wb.find("m:sheets", NS):
        tgt = rels[sh.get("{%s}id" % NS["r"])]
        tgt = tgt.lstrip("/") if tgt.startswith("/") else "xl/" + tgt
        rows = []
        for row in ET.fromstring(z.read(tgt)).iter("{%s}row" % NS["m"]):
            vals = {}
            for c in row.findall("m:c", NS):
                v = c.find("m:v", NS); t = c.get("t")
                if t == "inlineStr": x = "".join(e.text or "" for e in c.iter("{%s}t" % NS["m"]))
                elif v is None: x = ""
                elif t == "s": x = ss[int(v.text)]
                else: x = v.text
                vals[_col(c.get("r"))] = x
            rows.append([vals.get(i, "") for i in range(max(vals) + 1)] if vals else [])
        out[sh.get("name")] = rows
    return out
