"""Supplementary Figures 1-6 as one PDF: one page per figure, the figure (vector) with its legend below.
Inputs: the six figure PDFs (Supp Figure 0N - *.pdf) and "Supp Figures - legends.txt" in <supp dir>.
Needs: pdftocairo (poppler), Node.js with Playwright, qpdf.
Usage: python build_supp_figures_pdf.py <supp dir> <output.pdf>"""
import sys, os, re, glob, html, subprocess, tempfile, json

SUPP, OUT = sys.argv[1], os.path.abspath(sys.argv[2])
W = tempfile.mkdtemp()
leg = open(os.path.join(SUPP, "Supp Figures - legends.txt")).read()
blocks = {int(m.group(1)): m.group(0) for m in re.finditer(r"Supplementary Figure (\d)\.(?:.(?!\nSupplementary Figure \d\.))*", leg, re.S)}
LANDSCAPE = {2, 3}                                   # 11 x 8.5 in figures

def page(n, svg, land):
    text = " ".join(blocks[n].split())
    title, rest = text.split(". ", 1)
    width = "257mm" if land else "180mm"
    return f"""<!doctype html><html><head><meta charset="utf-8"><style>
@page{{size:A4 {'landscape' if land else 'portrait'};margin:12mm 15mm 12mm 15mm}}
body{{font-family:'Liberation Serif',Georgia,serif;font-size:9.5pt;line-height:1.4;color:#111;margin:0}}
.fig svg{{width:{width};height:auto;max-height:{'150mm' if land else '205mm'};display:block;margin:0 auto 4mm}}
.lg b{{font-family:'Liberation Sans',Arial,sans-serif}}
</style></head><body><div class="fig">{svg}</div><p class="lg"><b>{html.escape(title)}.</b> {html.escape(rest)}</p></body></html>"""

parts = []
for n in range(1, 7):
    pdf = glob.glob(os.path.join(SUPP, f"Supp Figure 0{n} - *.pdf"))[0]
    svgf = os.path.join(W, f"sf{n}.svg"); subprocess.run(["pdftocairo", "-svg", pdf, svgf], check=True)
    svg = open(svgf).read(); svg = svg[svg.index("<svg"):]
    h = os.path.join(W, f"p{n}.html"); open(h, "w").write(page(n, svg, n in LANDSCAPE)); parts.append((h, n in LANDSCAPE))
js = os.path.join(W, "r.js")
open(js, "w").write("""const { chromium } = require('playwright');
(async () => { const b = await chromium.launch(); const p = await b.newPage();
  for (const [h, land] of %s) { await p.goto('file://' + h); await p.pdf({ path: h.replace('.html', '.pdf'), preferCSSPageSize: true, printBackground: true }); }
  await b.close(); })();""" % json.dumps(parts))
env = dict(os.environ, NODE_PATH=subprocess.run(["npm", "root", "-g"], capture_output=True, text=True).stdout.strip())
subprocess.run(["node", js], check=True, env=env)
subprocess.run(["qpdf", "--empty", "--pages", *[h.replace(".html", ".pdf") for h, _ in parts], "--", OUT], check=True)
print("written", OUT)
