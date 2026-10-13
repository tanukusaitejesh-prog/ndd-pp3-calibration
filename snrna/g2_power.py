#!/usr/bin/env python3
"""Power and false-positive rate of the G2 tests (PREREG s7.4) on simulated genes. No real data are used.

Each simulated gene has 3 SNVs per position, per-site rates drawn from a lognormal with a Roulette-like spread,
and P(observed) = 1 - exp(-kappa * rate * selection), with kappa = 3 giving ~25-30% of SNVs observed (gnomAD v4.1
snRNA genes: 0.73 SNVs per nt, i.e. ~24%; disease-linked genes ~48%). "Misspecified" adds an unmodelled,
spatially autocorrelated rate error (AR(1) over positions, sd 0.5, lag-1 correlation 0.9) that the test's
expected counts do not know about, as a stand-in for Roulette missing local high-rate sites.

    python g2_power.py > results/g2_power_simulation.txt
"""
import numpy as np

import snrna


def simulate(L, region, factor, seed, misspecified=False, kappa=3.0):
    rng = np.random.default_rng(seed)
    n = np.repeat(np.arange(1, L + 1), 3)
    m = rng.lognormal(np.log(0.1), 0.8, len(n))
    true = m.copy()
    if misspecified:
        e = np.empty(L)
        e[0] = rng.normal(0, 0.5)
        for i in range(1, L):
            e[i] = 0.9 * e[i - 1] + np.sqrt(1 - 0.9 ** 2) * rng.normal(0, 0.5)
        true = m * np.exp(e[n - 1])
    sel = np.where(np.isin(n, region), factor, 1.0)
    obs = (rng.random(len(n)) < 1 - np.exp(-kappa * true * sel)).astype(float)
    p = 1 - np.exp(-snrna.fit_kappa(obs, m) * m)
    return n, obs, p


def rates(L, span, factor, sims, misspecified, seed0):
    region = list(range(span[0], span[1] + 1))
    hit = np.zeros(4)
    for s in range(sims):
        n, obs, p = simulate(L, region, factor, seed0 + s, misspecified)
        c = snrna.critical_region_test(n, obs, p, region)
        h = snrna.hotspot_test(n, obs, p, span)
        hit += [c["recovered"],                                                    # exact test alone (s7.4a rule)
                h["recovered_frozen"],                                             # frozen s7.4b rule
                c["recovered"] and h["top_decile_window"],                          # exact + top-decile window
                c["recovered"] and h["top_decile_window"] and h["p_circular"] < 0.05]  # proposed hybrid
    return hit / sims


def main():
    print("G2 test operating characteristics on simulated genes (fraction of simulations called 'recovered')")
    print("selection factor 1.0 = no depletion (false-positive rate); smaller = stronger depletion (power)\n")
    hdr = f"{'scenario':44s} {'factor':>6s} {'sims':>5s} {'exact':>7s} {'frozen':>7s} {'ex+top':>7s} {'hybrid':>7s}"
    cases = [("RNU4-2 critical region (13 nt / 145)", 145, None, [1.0, 0.05, 0.3]),
             ("RNU5B-1-like hotspot (8 nt / 117)", 117, (37, 44), [1.0, 0.02, 0.2, 0.5]),
             ("RNU2-2-like hotspot (67 nt / 191)", 191, (1, 67), [1.0, 0.5, 0.7])]
    for mis in (False, True):
        print(("MISSPECIFIED mutation model" if mis else "Correct mutation model"))
        print(hdr)
        for label, L, span, facs in cases:
            for f in facs:
                sims = 400 if f == 1.0 else 250
                if span is None:  # critical region is two segments: test it as the region set
                    region = snrna.CRITICAL
                    hit = np.zeros(4)
                    for s in range(sims):
                        rng_seed = 123000 + s + (50000 if mis else 0)
                        n, obs, p = simulate(L, region, f, rng_seed, mis)
                        c = snrna.critical_region_test(n, obs, p, region)
                        hit += [c["recovered"], np.nan, np.nan, c["recovered"] and c["p_circular"] < 0.05]
                    r = hit / sims
                else:
                    r = rates(L, span, f, sims, mis, 200000 + (50000 if mis else 0) + int(f * 1000))
                print(f"{label:44s} {f:6.2f} {sims:5d} " + " ".join("    n/a" if np.isnan(x) else f"{x:7.3f}" for x in r))
        print()
    print("Columns: exact = one-sided Poisson-binomial test vs gene-median expectation, p < 0.01 (frozen s7.4a);")
    print("frozen = top-decile window + circular-shift p < 0.01 (frozen s7.4b); ex+top = exact + top-decile window;")
    print("hybrid = exact + top-decile window + circular-shift p < 0.05 (for the critical region: exact + circular p < 0.05).")


if __name__ == "__main__":
    main()
