"""Assemble the public third cohort (data only, no analysis): GRCh37 de novo SNVs with sample, study and phenotype.
Sources: Jin et al. 2017 Nat Genet Table S9 (CHD cases) and denovo-db v1.6.1 (non-SSC) studies that are independent of
Kaplanis 2020 and Zhou 2022 and report complete coding calls (synonymous de novos >= 0.15 per proband)."""
import hashlib
import pandas as pd

X = "/mnt/x/ndd_gate_data" if __import__("os").path.exists("/mnt/x") else "X:/ndd_gate_data"
JIN = f"{X}/jin2017/jin2017_supplementary_tables.xlsx"
DDB = f"{X}/denovodb/denovo-db.non-ssc-samples.variants.v1.6.1.tsv.gz"
DDB_STUDIES = {"Fromer2014": "schizophrenia", "Gulsuner2013": "schizophrenia", "McCarthy2014": "schizophrenia",
               "Willsey2017": "tourette_syndrome", "epi4k2013": "epilepsy", "Rauch2012": "intellectualDisability"}
PROBANDS = {"Jin2017": 2645, "Fromer2014": 623, "Gulsuner2013": 105, "McCarthy2014": 57, "Willsey2017": 511,
            "epi4k2013": 264, "Rauch2012": 51}

j = pd.read_excel(JIN, sheet_name="S9", header=1)
j = j[(j["REF"].astype(str).str.len() == 1) & (j["ALT"].astype(str).str.len() == 1)]
j = pd.DataFrame({"study": "Jin2017", "phenotype": "congenital_heart_disease", "ndd_flag": j["NDD"].astype(str),
                  "sample": "Jin2017:" + j["Blinded ID"].astype(str), "chrom": j["CHROM"].astype(str),
                  "pos37": j["POS"].astype(int), "ref": j["REF"], "alt": j["ALT"]})

d = pd.read_csv(DDB, sep="\t", skiprows=1, low_memory=False).rename(columns={"#SampleID": "SampleID"})
d = d[d["StudyName"].isin(DDB_STUDIES) & (d["PrimaryPhenotype"] == d["StudyName"].map(DDB_STUDIES))]
d = d.drop_duplicates(["SampleID", "Chr", "Position", "Variant"])
ra = d["Variant"].astype(str).str.split(">", expand=True)
d = d.assign(ref=ra[0], alt=ra[1])
d = d[(d["ref"].str.len() == 1) & (d["alt"].str.len() == 1)]
d = pd.DataFrame({"study": d["StudyName"], "phenotype": d["PrimaryPhenotype"], "ndd_flag": "NA",
                  "sample": d["StudyName"] + ":" + d["SampleID"].astype(str), "chrom": d["Chr"].astype(str),
                  "pos37": d["Position"].astype(int), "ref": d["ref"], "alt": d["alt"]})
out = pd.concat([j, d], ignore_index=True)
out = out[out["chrom"].isin([str(i) for i in range(1, 23)])]
out["n_probands_study"] = out["study"].map(PROBANDS)
dst = f"{X}/public3/public3_dnv_snv_grch37.tsv"
__import__("os").makedirs(f"{X}/public3", exist_ok=True)
out.to_csv(dst, sep="\t", index=False)
print(out.groupby(["study", "phenotype"]).agg(carriers=("sample", "nunique"), snvs=("sample", "size"),
                                             probands=("n_probands_study", "first")).to_string())
print("total probands", sum(PROBANDS.values()), "SNVs", len(out), "md5", hashlib.md5(open(dst, "rb").read()).hexdigest())
