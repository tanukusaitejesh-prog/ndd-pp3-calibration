# A19 deviations log

1. **H1 baseline (found after running the pre-registered specification, 8 Oct 2026).** The pre-registered baseline
   (a = truth rate among PM1-applied, PP3-ABSENT SCVs) makes every rule an increment over "PP3 absent", so the learned
   credit includes the evidence carried by the absence of PP3 (learned Supporting credit 1.63 points). ACMG/AMP gives an
   absent code no points, and v13's credit estimand is ln LR(PP3 | PM1) relative to all PM1-applied SCVs. We report the
   pre-registered specification as primary and the ACMG-convention baseline (`make_rules_cv.py --acmg`) as a labelled
   sensitivity analysis. The H1a decision is the same under both (learned credit not better than the cap).

# A21 deviations
1. **Sensitivity analyses added after the pre-registered A21 run (8 Oct 2026).** Two problems surfaced: (a) one TP53
   score set (urn:mavedb:00001235-a-1, a dominant-negative-type readout) had inverted truth under the pre-registered
   nonsense/synonymous midpoint rule (LR(PP3) 0.33 < 1); (b) three F9 score sets had degenerate R (> 100) because a joint
   cell had no non-damaging variants (no smoothing was pre-specified). The pre-registered primary estimate is reported as
   primary; analyses excluding (a) and (a)+(b) are labelled sensitivity analyses.
