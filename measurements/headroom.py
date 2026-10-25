"""Two measurements on tree-sequence storage, framed as version-control questions.

1. Delta-base selection. Git packfiles search a window of candidate objects for the best
   delta base. A tree sequence always diffs a tree against its immediate predecessor.
   Is some earlier tree in a window a better base?

2. Content addressing. Git identifies objects by hashing their content, so identical content
   is stored once. tskit stores an edge once per contiguous span, but a parent-child
   relationship that disappears and later reappears costs a second record. How much
   non-contiguous recurrence is there?

Run:  python measurements.py
"""
import collections
import numpy as np
import msprime
import tskit

WINDOW = 50


def simulate(seed, n=200, L=5e6):
    return msprime.sim_ancestry(
        samples=n // 2, ploidy=2, sequence_length=L,
        recombination_rate=1e-8, population_size=10_000, random_seed=seed)


def content_addressing(ts):
    counts = collections.Counter((e.parent, e.child) for e in ts.edges())
    stored = sum(counts.values())
    distinct = len(counts)
    naive = sum(t.num_edges for t in ts.trees())
    return dict(distinct=distinct, stored=stored, naive=naive,
                saving=100 * (stored - distinct) / stored,
                repeated=100 * sum(1 for v in counts.values() if v > 1) / distinct,
                worst=max(counts.values()),
                compression=naive / stored)


def delta_base(ts, window=WINDOW):
    sets = [frozenset((t.parent(u), u) for u in t.nodes() if t.parent(u) != tskit.NULL)
            for t in ts.trees()]
    adj, best, offs = [], [], []
    for i in range(1, len(sets)):
        adj.append(len(sets[i] ^ sets[i - 1]))
        cand = [(len(sets[i] ^ sets[j]), i - j) for j in range(max(0, i - window), i)]
        b, off = min(cand)
        best.append(b); offs.append(off)
    adj, best, offs = np.array(adj), np.array(best), np.array(offs)
    return dict(adj=adj.mean(), best=best.mean(),
                reduction=100 * (1 - best.mean() / adj.mean()),
                nonadjacent=100 * (offs != 1).mean())


if __name__ == "__main__":
    for seed in (7, 23, 101):
        ts = simulate(seed)
        c, d = content_addressing(ts), delta_base(ts)
        print(f"\nseed {seed}: {ts.num_samples} samples, {ts.num_trees:,} trees, "
              f"{ts.num_edges:,} edges")
        print(f"  content addressing : {c['distinct']:,} distinct relationships stored as "
              f"{c['stored']:,} records")
        print(f"                       perfect dedup saves {c['saving']:.1f}%; "
              f"tskit already {c['compression']:.0f}x over naive ({c['naive']:,} records)")
        print(f"                       {c['repeated']:.1f}% of relationships recur, "
              f"worst recurs {c['worst']}x")
        print(f"  delta base         : mean edit vs previous {d['adj']:.2f}, "
              f"vs best in window({WINDOW}) {d['best']:.2f}")
        print(f"                       reduction {d['reduction']:.1f}%; "
              f"non-adjacent base better at {d['nonadjacent']:.1f}% of breakpoints")
