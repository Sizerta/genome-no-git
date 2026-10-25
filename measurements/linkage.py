"""Linkage disequilibrium is tree similarity - measured.

Two variants are correlated to the extent that their local genealogies agree:
  - identical carrier sets (the same lineage carries both) -> r^2 = 1
  - otherwise, correlation decays with the number of recombination breakpoints between
    them, because each breakpoint rewrites part of the tree

Pairs with identical carrier sets are reported separately rather than left inside the
distance bins: a lineage can persist across many breakpoints, and leaving those pairs in
inflates the distant bins and makes the decay look non-monotonic.

Run:  python linkage.py
"""
import collections
import numpy as np
import msprime

SEED, N, L, MAXDIST = 11, 300, 2e6, 50_000
MIN_FREQ, MIN_PAIRS = 0.05, 30


def load():
    ts = msprime.sim_ancestry(samples=N // 2, ploidy=2, sequence_length=L,
                              recombination_rate=1e-8, population_size=10_000,
                              random_seed=SEED)
    return msprime.sim_mutations(ts, rate=1e-8, random_seed=SEED + 1)


def main():
    ts = load()
    recurrent = {s.position for s in ts.sites() if len(s.mutations) > 1}
    G = ts.genotype_matrix()
    pos = np.array([s.position for s in ts.sites()])
    # common variants only: r^2 between rare variants is capped by frequency, not linkage
    freq = G.mean(1)
    keep = (freq > MIN_FREQ) & (freq < 1 - MIN_FREQ)
    G, pos = G[keep], pos[keep]
    tree_idx = []
    t = ts.first()
    for p in pos:
        while not (t.interval.left <= p < t.interval.right):
            t.next()
        tree_idx.append(t.index)
    tree_idx = np.array(tree_idx)
    print(f"{ts.num_samples} samples, {ts.num_trees:,} trees, {ts.num_sites:,} sites "
          f"({len(recurrent)} with >1 mutation), {len(pos):,} common variants")

    Gf = G.astype(np.float64)
    Z = (Gf - Gf.mean(1, keepdims=True)) / Gf.std(1, keepdims=True).clip(1e-12)
    bins, identical, imperfect = collections.defaultdict(list), [], 0
    for i in range(len(pos)):
        hi = np.searchsorted(pos, pos[i] + MAXDIST)
        for j in range(i + 1, hi):
            r2 = float((Z[i] @ Z[j]) / Z.shape[1]) ** 2
            if np.array_equal(G[i], G[j]):
                identical.append(r2)
            else:
                bins[min(abs(tree_idx[j] - tree_idx[i]), 21)].append(r2)
                if tree_idx[i] == tree_idx[j] and r2 > 0.99:
                    imperfect += pos[i] in recurrent or pos[j] in recurrent

    print(f"\n  {'relationship':>34} {'mean r^2':>9} {'n pairs':>9}")
    print("  " + "-" * 55)
    print(f"  {'identical carrier set (one lineage)':>34} {np.mean(identical):>9.3f} "
          f"{len(identical):>9,}")
    for d in [0, 1, 2, 3, 5, 10, 20, 21]:
        if len(bins.get(d, [])) > MIN_PAIRS:
            lab = ("same local tree" if d == 0 else
                   "more than 20 breakpoints apart" if d == 21 else
                   f"{d} breakpoint{'s' if d > 1 else ''} apart")
            print(f"  {lab:>34} {np.mean(bins[d]):>9.3f} {len(bins[d]):>9,}")

    print(f"\n  {'physical distance':>34} {'mean r^2':>9} {'n pairs':>9}")
    print("  " + "-" * 55)
    pb = collections.defaultdict(list)
    for i in range(len(pos)):
        hi = np.searchsorted(pos, pos[i] + MAXDIST)
        for j in range(i + 1, hi):
            pb[int((pos[j] - pos[i]) // 10_000)].append(
                float((Z[i] @ Z[j]) / Z.shape[1]) ** 2)
    for k in sorted(pb):
        if len(pb[k]) > MIN_PAIRS:
            print(f"  {f'{k*10}-{(k+1)*10} kb':>34} {np.mean(pb[k]):>9.3f} {len(pb[k]):>9,}")
    print("\n  Carrier sets that are identical give r^2 = 1 exactly. Where they do not,")
    print("  the cause is recurrent or back mutation at one of the two sites.")


if __name__ == "__main__":
    main()
