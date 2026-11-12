"""Python port of the dataviz skill's validate_palette.js (checks 2-5; same thresholds, Machado 2009 CVD sim)."""
import sys, math, itertools
BAND = {"light": (0.43, 0.77), "dark": (0.48, 0.67)}; CHROMA = 0.10; CVD_T, CVD_F = 8.0, 6.0; NORMAL = 15.0; CON = 3.0
SURF = {"light": "#fcfcfb", "dark": "#1a1a19"}
M = {"protan": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
     "deutan": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
     "tritan": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]]}
srgb = lambda h: [int(h.strip().lstrip("#")[i:i + 2], 16) / 255 for i in (0, 2, 4)]
s2l = lambda c: c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
lin = lambda h: [s2l(c) for c in srgb(h)]
def lum(h): r, g, b = lin(h); return 0.2126 * r + 0.7152 * g + 0.0722 * b
def contrast(a, b): x, y = sorted([lum(a), lum(b)], reverse=True); return (x + 0.05) / (y + 0.05)
def oklab_lin(c):
    r, g, b = c; cb = lambda v: math.copysign(abs(v) ** (1 / 3), v)
    l = cb(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b); m = cb(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b)
    s = cb(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b)
    return [0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s, 1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s]
def sim(h, k): r = lin(h); return [max(0, min(1, sum(M[k][i][j] * r[j] for j in range(3)))) for i in range(3)]
def dE(a, b, k=None):
    A = oklab_lin(sim(a, k) if k else lin(a)); B = oklab_lin(sim(b, k) if k else lin(b)); return 100 * math.dist(A, B)
def validate(p, mode="light", pairs="adjacent"):
    lo, hi = BAND[mode]; ok = True; out = []
    L = {c: oklab_lin(lin(c))[0] for c in p}; C = {c: math.hypot(*oklab_lin(lin(c))[1:]) for c in p}
    off = [(c, round(L[c], 3)) for c in p if not lo <= L[c] <= hi]; ok &= not off; out.append(("Lightness band", not off, off))
    lowc = [(c, round(C[c], 3)) for c in p if C[c] < CHROMA]; ok &= not lowc; out.append(("Chroma floor", not lowc, lowc))
    pl = list(itertools.combinations(range(len(p)), 2)) if pairs == "all" else [(i, i + 1) for i in range(len(p) - 1)]
    w = min((dE(p[i], p[j], k), k, p[i], p[j]) for k in ("protan", "deutan") for i, j in pl)
    st = "pass" if w[0] >= CVD_T else "floor" if w[0] >= CVD_F else "fail"; ok &= st != "fail"; out.append(("CVD separation", st, f"worst {w[2]}<->{w[3]} {w[0]:.1f} ({w[1]})"))
    n = min((dE(p[i], p[j]), p[i], p[j]) for i, j in pl); s2 = n[0] >= NORMAL; ok &= s2; out.append(("Normal-vision floor", s2, f"worst {n[1]}<->{n[2]} {n[0]:.1f}"))
    low = [(c, round(contrast(c, SURF[mode]), 2)) for c in p if contrast(c, SURF[mode]) < CON]; out.append(("Contrast vs surface", "relief" if low else "pass", low))
    return ok, out
if __name__ == "__main__":
    p = [c for c in sys.argv[1].split(",") if c]; mode = sys.argv[2] if len(sys.argv) > 2 else "light"; pairs = sys.argv[3] if len(sys.argv) > 3 else "adjacent"
    ok, out = validate(p, mode, pairs)
    for r in out: print(f"  {r[0]:22s} {r[1]!s:6s} {r[2]}")
    print("RESULT:", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
