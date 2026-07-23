#!/usr/bin/env python3
"""checkfig.py -- geometric sanity check for a compiled standalone figure.

Substitute for eyeballing: extracts every word's bounding box from the PDF
and reports (a) the figure extent, (b) any pair of words whose boxes
overlap by more than a tolerance, (c) words that fall outside the page.

Usage:  python3 checkfig.py file.pdf [expected_width_cm]
"""
import subprocess, sys, re, xml.etree.ElementTree as ET

PT_PER_CM = 72.0 / 2.54


def words(pdf):
    out = subprocess.run(["pdftotext", "-bbox", pdf, "-"],
                         capture_output=True, text=True).stdout
    # pdftotext emits raw control bytes for glyphs it cannot map to
    # Unicode (common with math fonts); those are illegal in XML.
    out = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", "?", out)
    root = ET.fromstring(out)
    ns = {"h": "http://www.w3.org/1999/xhtml"}
    pages = root.findall(".//h:page", ns)
    res = []
    for p in pages:
        pw = float(p.get("width")); ph = float(p.get("height"))
        for w in p.findall(".//h:word", ns):
            res.append((float(w.get("xMin")), float(w.get("yMin")),
                        float(w.get("xMax")), float(w.get("yMax")),
                        (w.text or "").strip()))
        return res, pw, ph
    return res, 0, 0


def overlap(a, b):
    ix = min(a[2], b[2]) - max(a[0], b[0])
    iy = min(a[3], b[3]) - max(a[1], b[1])
    if ix <= 0 or iy <= 0:
        return 0.0
    inter = ix * iy
    amin = min((a[2]-a[0])*(a[3]-a[1]), (b[2]-b[0])*(b[3]-b[1]))
    return inter / amin if amin > 0 else 0.0


def main():
    pdf = sys.argv[1]
    expected_cm = float(sys.argv[2]) if len(sys.argv) > 2 else None
    ws, pw, ph = words(pdf)
    print(f"  page: {pw/PT_PER_CM:.2f} cm x {ph/PT_PER_CM:.2f} cm, "
          f"{len(ws)} text items")
    if expected_cm:
        if pw/PT_PER_CM > expected_cm:
            print(f"  !! WIDER than the {expected_cm:.2f} cm target "
                  f"(will overflow the text block)")
        else:
            print(f"  OK width fits {expected_cm:.2f} cm target "
                  f"({100*pw/PT_PER_CM/expected_cm:.0f}% of it)")
    bad = []
    for i in range(len(ws)):
        for j in range(i+1, len(ws)):
            f = overlap(ws[i], ws[j])
            if f > 0.30:
                bad.append((f, ws[i][4], ws[j][4]))
    if bad:
        bad.sort(reverse=True)
        print(f"  !! {len(bad)} overlapping text pair(s):")
        for f, a, b in bad[:10]:
            print(f"       {f*100:4.0f}%  '{a}'  vs  '{b}'")
    else:
        print("  OK no overlapping text")
    off = [w for w in ws if w[0] < -1 or w[2] > pw+1 or w[1] < -1 or w[3] > ph+1]
    if off:
        print(f"  !! {len(off)} item(s) outside the page: "
              f"{[w[4] for w in off][:6]}")
    else:
        print("  OK all text inside the page")


if __name__ == "__main__":
    main()
