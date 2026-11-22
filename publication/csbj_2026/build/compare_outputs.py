"""Compare regenerated supplement files with the delivered ones (text content only; HTML tags, styles and whitespace ignored).
Usage: python compare_outputs.py <regenerated dir> <delivered dir> [file ...]
Prints IDENTICAL or the differing text lines for every file present in both folders."""
import sys, os, re, html, difflib
def text(p):
    s = open(p, encoding="utf8").read()
    if p.endswith(".html"):
        s = re.sub(r"<(style|script)[^>]*>.*?</\1>", " ", s, flags=re.S)
        s = re.sub(r'src="data:[^"]+"', "", s)
        s = re.sub(r"<(br|/p|/li|/tr|/h[1-6]|/div|/pre)\b[^>]*>", "\n", s)
        s = html.unescape(re.sub(r"<[^>]+>", " ", s))
    return [l for l in (re.sub(r"[ \t ]+", " ", x).strip() for x in s.splitlines()) if l]
new, old = sys.argv[1], sys.argv[2]
names = sys.argv[3:] or sorted(f for f in os.listdir(new) if f.endswith((".html", ".txt")))
bad = 0
for n in names:
    a, b = os.path.join(new, n), os.path.join(old, n)
    if not os.path.exists(b): print(f"[NOT DELIVERED] {n}"); continue
    d = list(difflib.unified_diff(text(b), text(a), "delivered", "regenerated", n=0, lineterm=""))
    if d: bad += 1; print(f"[DIFFERS] {n}"); print("\n".join(d[2:60]))
    else: print(f"[IDENTICAL] {n}")
sys.exit(1 if bad else 0)
