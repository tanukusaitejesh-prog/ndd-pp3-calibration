#!/usr/bin/env python3
"""Static, self-contained variant lookup page (one HTML file, no server) from the release scores and v4 thresholds.

    python make_lookup.py [--release release/snrna-vep-v2-final-<date>/scores.tsv.gz] [--out release/lookup.html]

A user picks a gene and types a transcript HGVS (n.4G>A) or a GRCh38 ID (chr11-62841806-C-T) and gets the frozen v4
score, its within-gene percentile, the PP3/BP4 evidence under the gene's own post hoc calibration (or the pooled one),
CADD for comparison, and the caveats. Research use; evidence labels are one line of evidence, never a classification.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STRENGTHS = ("Strong", "Moderate", "Supporting")


def thresholds(path):
    t = pd.read_csv(path, sep="\t")
    t = t[(t.prior == 0.10) & (t.score == "v4")]
    out = {}
    for r in t.itertuples():
        d = {f"{s}_{k}": getattr(r, f"{s}_{k}") for s in ("PP3", "BP4") for k in STRENGTHS}
        out[r.unit] = {k: (None if pd.isna(v) else float(v)) for k, v in d.items()}
        out[r.unit]["auroc"], out[r.unit]["n_P"], out[r.unit]["n_B"] = float(r.auroc), int(r.n_P), int(r.n_B)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", default=None)
    ap.add_argument("--thresholds", default=os.path.join(HERE, "results", "v4_calibration_thresholds.tsv"))
    ap.add_argument("--out", default=os.path.join(HERE, "release", "lookup.html"))
    a = ap.parse_args()
    rel = a.release or sorted(glob.glob(os.path.join(HERE, "release", "snrna-vep-v2-final-*", "scores.tsv.gz")))[-1]
    sc = pd.read_csv(rel, sep="\t")
    sc = sc[sc.gene.notna() & sc.v4.notna()]
    named = sc.groupby("gene").filter(lambda x: not str(x.gene.iloc[0]).startswith(("U", "ENSG")))
    named = named.assign(pct=named.groupby("gene").v4.rank(pct=True))
    hgvs = "n." + named.n.astype(int).astype(str) + named.ref_t.str.replace("U", "T") + ">" + named.alt_t.str.replace("U", "T")
    data = {}
    for g, x in named.assign(h=hgvs).groupby("gene"):
        data[g] = [[r.ID, r.h, round(float(r.v4), 4), round(float(r.pct), 4),
                    None if pd.isna(r.cadd_phred) else round(float(r.cadd_phred), 2)] for r in x.itertuples()]
    th = thresholds(a.thresholds)
    html = PAGE.replace("__DATA__", json.dumps(data, separators=(",", ":"))).replace("__TH__", json.dumps(th)) \
               .replace("__RELEASE__", os.path.basename(os.path.dirname(rel)))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(html)
    print(f"wrote {a.out}: {sum(len(v) for v in data.values())} SNVs in {len(data)} named genes; "
          f"calibrated units: {', '.join(th)} ({os.path.getsize(a.out) / 1e6:.1f} MB)")


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>snRNA variant lookup</title>
<style>
:root{--bg:#fcfcfb;--panel:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--muted:#8a8984;--line:#e4e3df;--accent:#2a78d6;
--pp3:#b42318;--bp4:#1a7f4b;--none:#52514e}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#1a1a19;--panel:#232322;--ink:#ffffff;--ink2:#c3c2b7;
--muted:#8f8e86;--line:#3a3936;--accent:#3987e5;--pp3:#f07a6e;--bp4:#5cc58f;--none:#c3c2b7}}
:root[data-theme="dark"]{--bg:#1a1a19;--panel:#232322;--ink:#ffffff;--ink2:#c3c2b7;--muted:#8f8e86;--line:#3a3936;--accent:#3987e5;
--pp3:#f07a6e;--bp4:#5cc58f;--none:#c3c2b7}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:760px;margin:0 auto;padding:28px 16px 48px}h1{font-size:22px;margin:0 0 4px}p.sub{color:var(--ink2);margin:0 0 22px}
.row{display:flex;gap:10px;flex-wrap:wrap;margin-bottom:14px}select,input,button{font:inherit;padding:9px 11px;border:1px solid var(--line);
border-radius:8px;background:var(--panel);color:var(--ink)}input{flex:1;min-width:200px}button{background:var(--accent);color:#fff;border:0;cursor:pointer}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:18px;margin-top:10px}
.big{font-size:28px;font-weight:600}.lab{color:var(--ink2);font-size:13px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:14px;margin-top:12px}
.ev{display:inline-block;padding:3px 10px;border-radius:99px;border:1px solid currentColor;font-weight:600;font-size:13px}
.PP3{color:var(--pp3)}.BP4{color:var(--bp4)}.none{color:var(--none)}.note{color:var(--ink2);font-size:13px;margin-top:16px}
.bar{height:8px;background:var(--line);border-radius:4px;position:relative;margin-top:6px}.bar i{position:absolute;top:-3px;width:3px;height:14px;background:var(--accent)}
table{border-collapse:collapse;width:100%;font-size:13px;margin-top:8px}td,th{border-bottom:1px solid var(--line);padding:5px 6px;text-align:left}
</style></head><body><main>
<h1>snRNA variant lookup</h1>
<p class="sub">Frozen v4 scores with ClinGen-calibrated PP3/BP4 evidence for spliceosomal snRNA genes. Release __RELEASE__. Research use.</p>
<div class="row"><select id="gene"></select><input id="q" placeholder="n.4G>A  or  chr11-62841806-C-T" autocomplete="off">
<button id="go">Look up</button></div>
<div id="out"></div>
<div class="note">Evidence uses the gene's own calibration where one exists (RNU2-2, RNU4-2, RNU4ATAC), otherwise the pooled
calibration, at prior 0.10 (Pejaver et al. 2022). It is one line of evidence (PP3/BP4), never a classification. For RNU4-2, the
measured saturation genome editing score should be used where available. CADD is shown for comparison only; it is not calibrated
for these genes.</div>
</main>
<script>
const DATA=__DATA__, TH=__TH__, ORDER=["Strong","Moderate","Supporting"];
const sel=document.getElementById("gene"), q=document.getElementById("q"), out=document.getElementById("out");
const pref=["RNU2-2","RNU4-2","RNU4ATAC","RNU5B-1","RNU12","RNU6ATAC"];
Object.keys(DATA).sort((a,b)=>((pref.indexOf(a)+1)||99)-((pref.indexOf(b)+1)||99)||a.localeCompare(b)).forEach(g=>{
  const o=document.createElement("option");o.value=g;o.textContent=g;sel.appendChild(o);});
function evidence(s,t){if(!t)return null;for(const k of ORDER){const v=t["PP3_"+k];if(v!==null&&v!==undefined&&s>=v)return["PP3",k];}
  for(const k of ORDER){const v=t["BP4_"+k];if(v!==null&&v!==undefined&&s<=v)return["BP4",k];}return["none",""];}
function esc(x){return String(x).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));}
function show(){
  const g=sel.value, s=q.value.trim().replace(/\s+/g,"").replace(/U/g,"T"); if(!s){out.innerHTML="";return;}
  const rows=DATA[g]||[]; const hit=rows.find(r=>r[0].toLowerCase()===s.toLowerCase()||r[1].toLowerCase()===s.toLowerCase());
  if(!hit){const m=s.match(/^n\.(\d+)/i), pos=m?m[1]:null;
    const near=pos?rows.filter(r=>r[1].match(/^n\.(\d+)/)[1]===pos).slice(0,6):[];
    out.innerHTML=`<div class="card">No SNV <b>${esc(s)}</b> in ${esc(g)}. Use transcript numbering (n.&lt;pos&gt;&lt;ref&gt;&gt;&lt;alt&gt;) or chr-pos-ref-alt (GRCh38).`+
    (near.length?`<table><tr><th>Did you mean</th><th>GRCh38</th></tr>${near.map(r=>`<tr><td>${esc(r[1])}</td><td>${esc(r[0])}</td></tr>`).join("")}</table>`:"")+`</div>`;return;}
  const [id,h,v4,pct,cadd]=hit; const own=TH[g], t=own||TH["pooled (all genes)"]; const ev=evidence(v4,t);
  const evh=ev?(ev[0]==="none"?`<span class="ev none">No PP3/BP4 at this score</span>`:`<span class="ev ${ev[0]}">${ev[0]} ${ev[1]}</span>`):"<span class='ev none'>not calibrated</span>";
  out.innerHTML=`<div class="card"><div class="lab">${esc(g)} ${esc(h)} &middot; ${esc(id)}</div>
  <div class="grid"><div><div class="lab">Evidence (${own?"gene-specific":"pooled"} calibration)</div><div style="margin-top:6px">${evh}</div></div>
  <div><div class="lab">v4 score</div><div class="big">${v4.toFixed(2)}</div></div>
  <div><div class="lab">Percentile within ${esc(g)}</div><div class="big">${Math.round(pct*100)}</div><div class="bar"><i style="left:calc(${(pct*100).toFixed(1)}% - 1px)"></i></div></div>
  <div><div class="lab">CADD v1.7 (reference only)</div><div class="big">${cadd===null?"&ndash;":cadd.toFixed(1)}</div></div></div>
  ${t?`<div class="note">Calibration unit: ${own?esc(g):"pooled (all genes)"}; ${t.n_P} P/LP vs ${t.n_B} population SNVs, AUROC ${t.auroc.toFixed(2)}.
  PP3 Supporting: ${t.PP3_Supporting===null?"not reached":"v4 &ge; "+t.PP3_Supporting.toFixed(2)}; BP4 Supporting: ${t.BP4_Supporting===null?"not reached":"v4 &le; "+t.BP4_Supporting.toFixed(2)}.</div>`:""}</div>`;}
document.getElementById("go").onclick=show; q.addEventListener("keydown",e=>{if(e.key==="Enter")show();}); sel.onchange=show;
</script></body></html>
"""

if __name__ == "__main__":
    main()
