"""A12.5: ClinGen Evidence Repository classifications recomputed under (a) uncapped, (b) PP3+PM1 capped at 4,
(c) conditional PP3 points (A10.1 tier-specific LR(PP3 | PM1-like), AlphaMissense DD, 2.406 scale, floored at 0)."""
import math, re, sys
import pandas as pd

src = sys.argv[1]
d = pd.read_csv(src, sep="\t", dtype=str).fillna("")
BASE = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BA": -8, "BS": -4, "BP": -1}
MOD = {"verystrong": 8, "strong": 4, "moderate": 2, "supporting": 1}
COND = {1: max(0.0, math.log(0.4469) / math.log(2.406)), 2: max(0.0, math.log(0.8994) / math.log(2.406)),
        4: math.log(2.634) / math.log(2.406)}


def pts(code):
    m = re.match(r"(PVS|PS|PM|PP|BA|BS|BP)(\d+)(?:_(\w+))?", code.strip())
    if not m:
        return None, 0
    v = BASE[m.group(1)]
    if m.group(3) and m.group(3).lower().replace("_", "") in MOD:
        v = MOD[m.group(3).lower().replace("_", "")] * (1 if v > 0 else -1)
    return m.group(1) + m.group(2), v


def cls(x):
    return "P" if x >= 10 else "LP" if x >= 6 else "VUS" if x >= 0 else "LB" if x >= -6 else "B"


rows = []
for r in d.itertuples():
    ev = {}
    for c in re.split(r"[,;]\s*", r[10]):
        k, v = pts(c)
        if k:
            ev[k] = v
    if "PP3" not in ev or "PM1" not in ev:
        continue
    other = sum(v for k, v in ev.items() if k not in ("PP3", "PM1"))
    pp3, pm1 = ev["PP3"], ev["PM1"]
    a = other + pp3 + pm1
    b = other + min(pp3 + pm1, 4)
    c = other + pm1 + COND[pp3 if pp3 in (1, 2) else 4] if pp3 > 0 else a
    rows.append({"gene": r[5], "assertion": r[9], "pp3": pp3, "pm1": pm1, "uncapped": cls(a), "capped": cls(b), "conditional": cls(c)})
x = pd.DataFrame(rows)
out = [f"A12.5: ClinGen Evidence Repository ({len(d)} classifications); PP3 and PM1 both applied: {len(x)}"]
for rule in ("capped", "conditional"):
    ch = x[x["uncapped"] != x[rule]]
    out.append(f"uncapped -> {rule}: {len(ch)} change class; " +
               ", ".join(f"{k} {v}" for k, v in (ch["uncapped"] + "->" + ch[rule]).value_counts().items()))
out.append("published assertion counts: " + ", ".join(f"{k} {v}" for k, v in x["assertion"].value_counts().items()))
out.append("PP3 strength among these: " + ", ".join(f"{k} pt {v}" for k, v in x["pp3"].value_counts().sort_index().items()))
out.append("top genes: " + ", ".join(f"{k} {v}" for k, v in x["gene"].value_counts().head(8).items()))
open("results_fullstudy_a12_5_summary.txt", "w").write("\n".join(out) + "\n")
x.to_csv("a12_5_erepo_pp3_pm1.tsv", sep="\t", index=False)
print("\n".join(out))
