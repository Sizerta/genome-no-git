"""Simulate a small, legible ARG with msprime and export it as JSON for the explainer.
Fixed leaf order across all local trees so that topology changes are readable when the
genome position moves: samples stay put, branches move."""
import json, sys
import numpy as np
import msprime, tskit

NE, L, N = 10_000, 50_000, 8

def layout(tree, ts, leaf_x):
    """x: leaves fixed, internal = mean of children. y: node time."""
    xs, ys = {}, {}
    for u in tree.nodes(order="postorder"):
        ys[u] = ts.node(u).time
        ch = tree.children(u)
        xs[u] = leaf_x[u] if not ch else float(np.mean([xs[c] for c in ch]))
    return xs, ys

def build(seed):
    ts = msprime.sim_ancestry(samples=N // 2, ploidy=2, sequence_length=L,
                              recombination_rate=2.5e-9,
                              population_size=NE, random_seed=seed)
    ts = msprime.sim_mutations(ts, rate=5.8e-9, random_seed=seed + 1)
    return ts

best = None
for seed in range(1, 400):
    ts = build(seed)
    nt, ns = ts.num_trees, ts.num_sites
    if 5 <= nt <= 7 and 24 <= ns <= 40:
        # prefer examples where topologies genuinely differ between adjacent trees
        tops = [t.rank() for t in ts.trees()]
        changes = sum(1 for a, b in zip(tops, tops[1:]) if a != b)
        score = changes
        if best is None or score > best[0]:
            best = (score, seed, nt, ns)
print("chosen:", best, file=sys.stderr)
score, seed, nt, ns = best
ts = build(seed)

leaf_x = {u: i for i, u in enumerate(ts.samples())}
trees = []
prev_parent = None
for t in ts.trees():
    xs, ys = layout(t, ts, leaf_x)
    nodes = sorted(xs)
    parent = {int(u): (int(t.parent(u)) if t.parent(u) != tskit.NULL else None) for u in nodes}
    muts = []
    for site in t.sites():
        for m in site.mutations:
            muts.append(dict(pos=float(site.position), node=int(m.node),
                             derived=m.derived_state, anc=site.ancestral_state))
    changed = (None if prev_parent is None else
               sorted(int(u) for u in set(parent) | set(prev_parent)
                      if parent.get(u) != prev_parent.get(u)))
    trees.append(dict(
        left=float(t.interval.left), right=float(t.interval.right),
        nodes=[dict(id=int(u), x=xs[u], y=ys[u], is_sample=bool(t.is_sample(u)),
                    parent=parent[int(u)], n_below=int(t.num_samples(u))) for u in nodes],
        mutations=muts, tmrca=float(max(ys.values())),
        total_branch_length=float(t.total_branch_length),
        changed_nodes=changed))
    prev_parent = parent

G = ts.genotype_matrix()
sites = [dict(pos=float(s.position), idx=i) for i, s in enumerate(ts.sites())]
# which node does each site's mutation sit on, and which samples carry the derived allele
site_info = []
for i, s in enumerate(ts.sites()):
    m = s.mutations[0]
    tr = next(j for j, t in enumerate(trees) if t["left"] <= s.position < t["right"])
    carriers = [int(k) for k in range(N) if G[i, k] == 1]
    site_info.append(dict(idx=i, pos=float(s.position), node=int(m.node),
                          tree=tr, carriers=carriers, freq=len(carriers) / N))

out = dict(n_samples=N, seq_len=L, Ne=NE, seed=seed,
           sample_ids=[int(u) for u in ts.samples()],
           breakpoints=[float(b) for b in ts.breakpoints()],
           trees=trees, sites=site_info,
           genotypes=[[int(v) for v in row] for row in G],
           num_trees=int(ts.num_trees), num_sites=int(ts.num_sites),
           num_edges=int(ts.num_edges), num_nodes=int(ts.num_nodes))
json.dump(out, open("/home/claude/argviz/arg.json", "w"))
print(f"seed={seed} trees={nt} sites={ns} nodes={ts.num_nodes} edges={ts.num_edges}")
print("breakpoints:", [round(b) for b in ts.breakpoints()])
for i, t in enumerate(trees):
    print(f"  tree {i}: [{t['left']:.0f},{t['right']:.0f})  tmrca={t['tmrca']:.0f} "
          f"muts={len(t['mutations'])} changed_nodes={t['changed_nodes']}")
