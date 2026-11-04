# Verified timeline of freezes and results (Supplementary Table S1)

Verified on 4 November 2026. Bitcoin times are the timestamps of the earliest block attesting each manifest's OpenTimestamps proof. Proofs were upgraded with `ots upgrade` on 4–4 November 2026, and block times were read from blockstream.info. Git times are commit times in UTC. A Bitcoin anchor proves the manifest existed by that time. A git time alone is weaker evidence: commit times are set by the committer, and independent confirmation rests on GitHub's push record.

| Event | Evidence | Time (UTC) |
|---|---|---|
| Pre-registration (`PREREG.md`) frozen | `prereg_manifest.txt.ots`, Bitcoin block 971237 | 2026-10-13 11:59 |
| Analysis code frozen | `code_manifest.txt.ots`, block 971241 | 2026-10-13 12:20 |
| Amendment A1 (G2 hotspot rule, simulation only) | block 971253 | 2026-10-13 14:07 |
| G2 mutation model frozen | `g2_freeze_manifest.txt.ots`, block 971255 | 2026-10-13 14:11 |
| Amendment A2 (v2 model specification) and A2b | block 971283 | 2026-10-13 18:33 |
| Amendment A4 (analysis plan for open endpoints) | block 971354 | 2026-10-14 05:07 |
| v4 code committed | git f4bd48a | 2026-10-14 05:49 |
| **A5: v4 code and DGX output hashes frozen** | `amendment_A5_manifest.txt.ots`, **block 971362** | **2026-10-14 05:58** |
| **A6: inheritance rule frozen** | `amendment_A6_manifest.txt.ots`, **block 971363** | **2026-10-14 06:04** |
| Blind (sealed) results committed | git 4081c92 | 2026-10-14 06:44 |
| A7 frozen (code) → scores frozen → evaluated | git 93682d0 → 55033c3 → 1d4086c | 09:47 → 10:10 → 10:20 |
| A8 frozen → labels committed → evaluated | git 4fdf25b → e1bb6e4 → 68710b0 | 10:23 → 14:13 → 14:16 |
| A9 frozen (positives) → evaluated | git 5b2145e → 4de7c80 | 10:32 → 14:22 |
| A10 frozen → missense evaluated | git a411288 → 13608ab | 15:16 → 17:15 |
| A11 frozen → evaluated | git e61a8f6 → 0de2637 | 10-14 15:18 → 10-15 06:13 |
| A12/A13 frozen → A13 evaluated | git cc399a8 → b9cec28 | 10-15 04:29 → 17:50 |
| A14/A15 frozen → A15 simulation run | git 35fdd9d → 4096c06 | 10-15 06:28 → 08:10 |
| A16 frozen (not yet run) | git b2818c8 | 10-15 18:07 |

## Notes
- **Sealing and decryption** of the RNU2-2 and RNU5B-1 labels are logged on the A100 machine in `sealed/SEALED_MANIFEST.jsonl`. That file must be added to the release before submission. It holds event times and hashes, not labels or the passphrase.
- **The A5 manifest records the DGX output hashes** (`v4_scores.tsv.gz`, `v4_model.json`) as 62 of 64 hex characters, transcribed from a screenshot. The full hashes are in `results_v4/v4_manifest.txt` on the DGX, and `v4_blind_eval.py` checks them before scoring. The 62-character suffix is still a unique fingerprint (248 bits).
- **The code hashes in the A5 and A6 manifests still match** the committed `v4.py`, `v4_blind_eval.py`, `inheritance_blind_eval.py` and `final_eval.py` (checked 4 November 2026).
- **Amendments A7–A16 have no OpenTimestamps proof in the repository**, so their order rests on git and GitHub push history.
