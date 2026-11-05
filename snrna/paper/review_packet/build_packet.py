#!/usr/bin/env python3
"""Assemble the external review packets (internal rating removed).   python paper/review_packet/build_packet.py"""
import os, re, zipfile
HERE = os.path.dirname(os.path.abspath(__file__)); P = os.path.dirname(HERE); ROOT = os.path.dirname(P)
rd = lambda f: open(f).read().strip()
ver = rd(os.path.join(P, "VERIFICATION_2026-11-04.md"))
ver = ver[:ver.index("## 7. Rating")].rstrip() + "\n\n## 7. Rating\n[Removed from the review packet so that the external rating is independent.]\n"
docs = [("MANUSCRIPT", rd(os.path.join(P, "MANUSCRIPT.md"))), ("SUPPLEMENT", rd(os.path.join(P, "SUPPLEMENT.md"))),
        ("TIMELINE", rd(os.path.join(P, "TIMELINE.md"))), ("VERIFICATION REPORT (rating removed)", ver),
        ("REFERENCES", rd(os.path.join(P, "REFERENCES.md")))]
prereg = ("PRE-REGISTRATION (PREREG.md: frozen plan, deviation log §14, amendments §15)", rd(os.path.join(ROOT, "PREREG.md")))
prompt = rd(os.path.join(HERE, "REVIEW_PROMPT.md"))
def build(name, items, note):
    out = [prompt, "\n\n---\n\n# MATERIALS UNDER REVIEW\n", note]
    for n, (title, body) in enumerate(items, 1):
        body = re.sub(r"^(#{1,5}) ", lambda m: "#" + m.group(1) + " ", body, flags=re.M)
        out.append(f"\n\n---\n\n# DOCUMENT {n}: {title}\n\n" + body)
    out.append("\n\n---\n\n# END OF MATERIALS\nNow complete tasks 1–8 from the instructions at the top.\n")
    t = "".join(out); open(os.path.join(HERE, name), "w").write(t); return len(t.split())
w1 = build("REVIEW_PACKET_FULL.md", docs + [prereg], "Six documents follow: manuscript, supplement, timeline, verification report, references and the full pre-registration.")
w2 = build("REVIEW_PACKET_COMPACT.md", docs, "Five documents follow: manuscript, supplement, timeline, verification report and references. The full pre-registration and deviation log are omitted for length; ask for them if you need a specific entry.")
with zipfile.ZipFile(os.path.join(HERE, "review_packet_files.zip"), "w", zipfile.ZIP_DEFLATED) as z:
    z.writestr("files/00_REVIEW_PROMPT.md", prompt + "\n")
    z.writestr("files/VERIFICATION_REPORT_no_rating.md", ver)
    for f in ("MANUSCRIPT.md", "SUPPLEMENT.md", "TIMELINE.md", "REFERENCES.md"):
        z.write(os.path.join(P, f), "files/" + f)
    z.write(os.path.join(ROOT, "PREREG.md"), "files/PREREG.md")
print(f"full {w1} words, compact {w2} words")
