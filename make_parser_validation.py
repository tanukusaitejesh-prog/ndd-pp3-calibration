"""Validation of the free-text ACMG/AMP code parser used for laboratory submissions (review 18). Ground truth: the
structured 'Applied Evidence Codes (Met)' of ClinGen Evidence Repository classifications. Test text: the free-text
'Summary of interpretation' of the same classifications. Reports, for codes mentioned in the summary, precision and
recall of the parser against the structured codes, overall and for PP3, PM1 and PP2 (presence and strength).
Run: python make_parser_validation.py work_erepo/erepo_2026-10-07.tsv"""
import re, sys
import pandas as pd

BASE = {"PVS": 8, "PS": 4, "PM": 2, "PP": 1, "BA": -8, "BS": -4, "BP": -1}
MOD = {"verystrong": 8, "strong": 4, "moderate": 2, "supporting": 1, "standalone": 8}
CODE = re.compile(r"\b(PVS|PS|PM|PP|BA|BS|BP)([1-9])(?:[_ ]?\(?(very[ _]?strong|strong|moderate|supporting|stand[ _-]?alone)\)?)?\b", re.I)
NEG_A = re.compile(r"^(?:/[A-Z]{2,3}\d)?\W*(?:criteria\s+|criterion\s+)?(?:was\s+|is\s+|were\s+|are\s+)?(?:not\b|n/a|unmet|absent\b)", re.I)
NEG_B = re.compile(r"(?:\b(?:not|no|without|neither|nor)\b|does\s+not\s+meet|do\s+not\s+meet|not\s+meet|insufficient\s+for|below\s+the\s+threshold\s+for)[^.;]{0,40}$", re.I)


def parse_free(txt):
    ev = {}
    for m in CODE.finditer(txt or ""):
        if NEG_A.match(txt[m.end():m.end() + 25]) or NEG_B.search(txt[max(0, m.start() - 60):m.start()]):
            continue
        fam, num, mod = m.group(1).upper(), m.group(2), m.group(3)
        v = BASE[fam]
        if mod:
            v = MOD[re.sub(r"[ _-]", "", mod.lower())] * (1 if v > 0 else -1)
        ev[fam + num] = v
    return ev


def parse_struct(txt):
    ev = {}
    for c in re.split(r"[,;]\s*", txt):
        m = re.match(r"(PVS|PS|PM|PP|BA|BS|BP)(\d+)(?:_(\w+))?", c.strip())
        if not m:
            continue
        v = BASE[m.group(1)]
        if m.group(3) and m.group(3).lower().replace("_", "") in MOD:
            v = MOD[m.group(3).lower().replace("_", "")] * (1 if v > 0 else -1)
        ev[m.group(1) + m.group(2)] = v
    return ev


d = pd.read_csv(sys.argv[1], sep="\t", dtype=str).fillna("")
d = d[d["Retracted"].str.lower() != "true"]
F = d["Summary of interpretation"].map(parse_free)
S = d["Applied Evidence Codes (Met)"].map(parse_struct)
NOTMET = d["Applied Evidence Codes (Not Met)"].map(parse_struct)
has = F.map(len) > 0
out = [f"Classifications {len(d)}; summaries mentioning >= 1 code {int(has.sum())}"]
tp = fp = fn = 0
fp_notmet = 0
per = {k: [0, 0, 0, 0, 0] for k in ("PP3", "PM1", "PP2")}       # tp, fp, fn, strength agree, strength compared
for f, s_, nm in zip(F[has], S[has], NOTMET[has]):
    for c in f:
        if c in s_:
            tp += 1
        else:
            fp += 1
            fp_notmet += c in nm
    fn += sum(1 for c in s_ if c not in f)
    for k in per:
        if k in f and k in s_:
            per[k][0] += 1; per[k][4] += 1; per[k][3] += f[k] == s_[k]
        elif k in f:
            per[k][1] += 1
        elif k in s_:
            per[k][2] += 1
out.append(f"All codes (summaries that mention codes): precision {tp / (tp + fp):.3f}, recall {tp / (tp + fn):.3f} "
           f"(TP {tp}, FP {fp} of which {fp_notmet} are listed as 'not met', FN {fn})")
for k, (a, b, c, agree, comp) in per.items():
    out.append(f"{k}: precision {a / max(a + b, 1):.3f}, recall {a / max(a + c, 1):.3f}, strength agreement {agree / max(comp, 1):.3f} "
               f"(TP {a}, FP {b}, FN {c})")
open("results_fullstudy_backup/parser_validation_summary.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
