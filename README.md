# The genome has no `.git` directory

**Population genetics for people who already think in commits.**

There is a field that asks a question you already know how to think about: given only the
files in the latest checkout, reconstruct the commit history that produced them. No
`.git` directory. No log. Just the working tree.

That field is population genetics, and the reconstructed history is called an **ancestral
recombination graph** — an ARG. If you write software, you have most of the intuition
already. It is in the wrong vocabulary, which is the only thing standing between you and
a literature that is currently one of the most interesting in genomics.

This piece maps the two worlds onto each other, uses the map to explain the quantities
the field actually cares about — linkage, relatedness, allele frequency — and then spends
the second half on the places the map fails. The failures are the useful part. Each one
is a concept you would otherwise have to learn cold.

There is an **[interactive version of the example](https://sizerta.github.io/genome-no-git/)** —
eight genomes, 50 kb, real simulated data. Open it in another tab. Everything below
refers to it.

---

## The mapping

| Git | Population genetics |
| --- | --- |
| a checked-out file | a sampled genome |
| commit | an ancestral genetic state |
| diff between commits | a mutation |
| `git blame` | which branch a variant arose on |
| branch | a lineage of inheritance |
| merge commit | recombination |
| `git log -- one/file` | the local tree for one region of the genome |
| the commit graph (a DAG) | the ARG |
| reconstructing history from a checkout | ARG inference |

The table is the cheap part. Anyone can generate it. What follows is the part that
earns its keep.

---

## One coordinate, one tree

Pick any single base in a genome. For a set of sampled genomes, there is exactly one
genealogy at that base: who shares an ancestor with whom, and how long ago. Draw it and
you get a tree, samples at the bottom, time running upward.

In the interactive example, eight haplotypes have a genealogy whose deepest node sits
around 43,000 generations back. Most pairs find a common ancestor far more recently than
that. The time axis defaults to a square root for a reason: switch it to linear and
almost everything collapses toward zero. That skew is not a drawing artefact. It is what
the coalescent process produces, and noticing it is most of the intuition for why
population genetics talks about deep time so much — the last few coalescences take
longer than everything before them combined.

## A mutation is an edge

Put a mutation on a branch. Every sample below that branch inherits it. Every sample
elsewhere in the tree does not.

So the set of genomes carrying a variant is exactly the set of leaves under one edge.
That carrier set is a column of a genotype matrix. If you have spent any time with
variant data, you have been looking at edges of trees without knowing it.

Two consequences fall out with no algebra:

- **Allele frequency is tree position.** A mutation high in the tree has many descendants
  and is common. A mutation on the last twig before a single sample is a singleton.
- **`git blame` has a precise analogue.** Asking which branch a variant arose on is the
  same question as asking which commit introduced a line.

Two variants from the simulated example, to make that concrete:

| | position | sits on | carried by | frequency |
| --- | --- | --- | --- | --- |
| a common variant | 1,194 bp | the branch above node 15 | haplotypes 0, 2, 3, 4, 5, 7 | 0.750 |
| a singleton | 2,940 bp | the branch above node 6 | haplotype 6 | 0.125 |

Nothing distinguishes these two mutations except where they landed. One fell on a branch
with six descendants, one on a twig with a single descendant. That is the entire
explanation of the allele frequency spectrum, which is otherwise one of the first things
you have to take on faith.

In the interactive version, hovering a column of the genotype matrix lights up the branch
it sits on and the leaves that carry it. That single gesture is the whole idea.

## Different files, different histories

Here is where the analogy stops being decorative and starts doing work.

`git log -- src/parser.c` gives you a different history from `git log -- README.md`. Same
repository, different paths, different sequences of commits touching them.

Genomes are the same, and for a concrete physical reason. Your chromosome was not copied
from one parent's chromosome. It was assembled by cutting and splicing two of them —
recombination — roughly once or twice per chromosome per generation. So the first stretch
of a chromosome has one ancestral history and the next stretch has another.

**There is no single tree for a genome.** There is a sequence of local trees, one per
non-recombining interval, called marginal trees. The interactive example has seven of
them across 50 kb. Drag the marker across a breakpoint and the topology rearranges.

A recombination event is a merge commit in the strict sense: this region's ancestry comes
from one parent lineage, that region's from another. The ARG is the graph that holds all
the local trees at once, and its defining feature is that some nodes have two parents.

## Why tree sequences compress

Drag across a breakpoint in the example and watch what moves. Most of the tree does not.
Seven local trees there are built from 20 nodes, and neighbouring trees differ at only
four or five of them.

So you do not store seven trees. You store the edits. This is the insight behind
[`tskit`](https://tskit.dev), the library that makes genealogies for hundreds of
thousands of samples tractable: an edge is recorded once, with the genomic interval over
which it exists, instead of once per tree it appears in.

The effect is not subtle. In a simulation of 200 sampled haplotypes over 5 Mb
(parameters and code in [`measurements/`](measurements/)):

```
trees: 10,746     nodes: 7,292     edge records: 39,208
storing every tree's edges separately: 4,276,908 records
```

That is a **109× reduction**, and it comes from nothing more than noticing that
recombination edits a tree instead of replacing it.

---

## Linkage is tree similarity

Two variants in a genome are usually correlated: knowing someone's allele at one tells
you something about the other. This is linkage disequilibrium, and it is the quantity
most of statistical genetics is built on. Every GWAS hit you have ever read about is a
marker standing in for something nearby it is correlated with.

The genealogical picture explains it in one sentence. **Two mutations on the same branch
are carried by exactly the same samples, so they are perfectly correlated.** Not
approximately. Identically, by construction, at any physical distance.

As the genealogy is rewritten by recombination, the two carrier sets drift apart and the
correlation decays.

That is a measurable claim, so here it is measured — 300 sampled haplotypes over 2 Mb,
common variants only, grouped by how many recombination breakpoints separate the two
sites ([`measurements/linkage.py`](measurements/linkage.py)):

```
relationship                          mean r²    n pairs
identical carrier set (one lineage)     1.000        930
same local tree                         0.256        937
1 breakpoint apart                      0.288      1,123
2 breakpoints apart                     0.301      1,011
5 breakpoints apart                     0.292      1,142
10 breakpoints apart                    0.214      1,081
20 breakpoints apart                    0.162      1,069
more than 20 breakpoints apart          0.090     91,151
```

Identical carrier sets give **r² = 1.000**, exactly, as predicted. Beyond that the first
few bins are indistinguishable from each other — a handful of breakpoints need not disturb
the part of the tree your two variants sit on — and the decay only becomes clear past
about five.

The familiar version of this plot uses physical distance, and it is the same phenomenon
seen through a blurrier lens:

```
0-10 kb    0.252        30-40 kb   0.075
10-20 kb   0.130        40-50 kb   0.060
20-30 kb   0.095
```

Physical distance is a proxy. Genealogical distance is the thing itself. Two sites 40 kb
apart with no breakpoint between them are perfectly linked; two sites 2 kb apart with a
breakpoint between them need not be. LD decay curves are averages over a variable that is
not being measured directly.

One detail worth keeping, because it is the analogy failing in public. The r² for
identical carrier sets came out at 0.996 before I looked closely. Three pairs out of
hundreds were not perfectly correlated, and all three involved a site carrying more than
one mutation — the same base mutating twice. That is homoplasy, break 5 below,
demonstrating itself in a measurement designed for something else.

## Relatedness is a local quantity

Pick two haplotypes and ask how long ago they shared an ancestor. In the interactive
example, haplotypes 2 and 4 coalesce **712 generations ago** in one region and **34,688
generations ago** in another. A 49-fold swing, across 50 kb, between the same two
sequences. Near-siblings at one locus, distant cousins at the next.

The genome-length-weighted average is 18,789 generations, and that single number is what
a relatedness matrix stores for the pair — the matrix you hand to a mixed model when you
fit heritability or run an association test with a random effect for kinship.

So the kinship coefficient is not a property of two individuals. It is a genome-wide
average over a quantity that varies enormously along the chromosome, and the averaging
throws away structure you can see directly in the panel at the bottom of the interactive
page.

This is not a pedantic point. A 2025 preprint on fitting linear mixed models directly to
ARGs observes that narrow-sense heritability is not a single number for diploids, because
self-relatedness differs between individuals. That claim sounds strange until you have
watched one pair's coalescence time jump by a factor of 49 across a short stretch of
sequence, and then it sounds obvious.

## What this buys you

Here is the reason to care, stated plainly.

Allele frequency, linkage disequilibrium, relatedness, identity by descent, population
structure, admixture, the signature of natural selection, imputation of ungenotyped
variants, and the heritability of a trait are usually taught as separate topics with
separate estimators. Every one of them is a summary statistic of the genealogy.

If you had the true ARG, you would not estimate linkage disequilibrium — you would read it
off. You would not estimate kinship — you would compute it. The reason these are separate
subfields with separate methods is historical: nobody could get at the genealogy, so the
field built a different estimator for each shadow the genealogy casts.

That is what has changed in the last few years, and why this literature is worth your
attention. Inferring ARGs at the scale of hundreds of thousands of genomes became
possible, which turns a collection of estimators into one object that generates all of
them.

---

## Where the analogy breaks

Everything above makes population genetics sound like a solved engineering problem. It is
not, and the reasons are exactly the places the git picture fails. Each one is a concept
worth having.

### 1. There is no `.git`

Git does not reconstruct history. It records it, at the moment the work happens, and
retrieves it later. Biology records nothing. Every ancestral genome is gone. You have the
working tree and a theory of how copying works.

This is the whole field in one sentence. Population genetics is not parsing a log, it is
statistical inference of a log that was never written. Every method downstream of this is
shaped by that fact.

### 2. There are no timestamps

A commit carries its date. A mutation does not. You infer when a branch existed from how
many mutations accumulated along it, because mutations arrive at a roughly steady rate —
the molecular clock. Branch lengths in all those trees are not measured. They are
estimated from counting.

### 3. You can only walk backward from `HEAD`

Git history can be read in either direction. Genealogy can only be read one way: from the
present, upward. You have the leaves and nothing else, so every model has to run backward
in time. That is why the central object is called the **coalescent** — a process defined
by lineages merging as you go back, rather than splitting as you go forward.

### 4. Branches with no surviving descendants are gone forever

A deleted git branch is usually still recoverable from the reflog, at least until it is
garbage collected. A lineage that left no descendants in your sample is unrecoverable,
permanently and in principle — there is no reflog and no grace period. You can only ever reconstruct the
history ancestral to the samples you hold. Everything else is invisible, and no amount of
sequencing the living recovers it.

### 5. The same diff can happen twice

Git identifies a change uniquely. Biology does not. The same base can mutate to the same
alternative independently in two lineages, or mutate and then mutate back.

This matters more than it sounds. A great deal of convenient reasoning assumes every
mutation is unique and never reverts — the infinite sites assumption. Where it fails you
get homoplasy, and two genomes can share a variant without sharing the ancestor you would
infer from it. Detecting those conflicts is old, central machinery in the field.

### 6. Many logs explain the same checkout

This is the one to sit with.

Given a set of genomes, there is not one ARG consistent with them. There is an enormous
space of them, and the data does not single one out. Recombination makes it worse:
the number of histories that produce the same observed sequences grows combinatorially.

Most published analyses infer one ARG and then use it as though it were the truth. A 2025
review in *Nature Reviews Genetics* names this directly: treating a single inferred ARG
as the true history ignores that uncertainty and can bias the results, and the field does
not yet have well-calibrated inference. This is the open problem, and it is the kind a
software person recognises instantly — a point estimate being passed around as a fact.

### 7. Nothing in git is more likely to survive because it is good

This is the largest omission in the whole analogy, and it is worth stating on its own.

Git has no notion that one commit is more likely to persist than another because of what
it contains. Every commit in the history is there because somebody wrote it, not because
it earned its place.

Biology's central mechanism is exactly the thing git lacks. Lineages carrying some
variants leave more descendants than lineages carrying others, *because of what they
carry*. That is natural selection, and it bends the genealogy: a beneficial variant
sweeping through a population drags its neighbours with it and leaves a region whose local
trees are shallow and star-shaped, quite unlike the trees a neutral process produces.

Everything above — the coalescent, the branch-length intuition, the clock — describes
neutral evolution, where variants are passengers. That assumption is the backdrop against
which selection is detected: you characterise what the genealogy should look like if
nothing were being selected, then look for regions where it does not.

Worth flagging because it propagates. The ARG-based mixed model work notes that under
selection the genetic component no longer reduces to an intercept, and conjectures that
the resulting structure could produce false discoveries. Methods built on neutral
assumptions inherit those assumptions wherever they are applied.

---

## Where git's solutions do not transfer

An analogy this close invites a question: can the engineering move across? Git has twenty
years of tuning in it. Some of that should be worth something.

I tested the two most obvious candidates. Both came back flat, and the reasons are
instructive.

**Delta base selection.** Git's packfiles do not diff an object against the previous one.
They search a window of candidates for the one that produces the smallest delta. Tree
sequences diff blindly against the adjacent tree. So: within a 50-tree window, is some
earlier tree a better base?

```
mean edit size vs previous tree        : 7.22
mean edit size vs best in 50-tree window: 7.20   (0.4% better)
breakpoints where a non-adjacent base wins: 0.7%
```

Essentially nothing. And the reason is a genuine insight about the data rather than a
disappointment. **A commit can touch any file anywhere, so git has to search. A
recombination event edits the tree exactly where you are standing.** The adjacent tree is
almost always already the nearest neighbour in edit distance. Locality is built into the
biology, so the search has nothing to find.

**Content addressing.** Git identifies objects by hashing their content, so identical
objects are stored once no matter where they appear. `tskit` stores an edge once per
contiguous interval, but an edge that disappears and later reappears becomes two records.
Git would store that once.

```
distinct (parent, child) relationships : 34,476
actual edge records stored             : 39,208
theoretical saving                     : 12.1%
```

Twelve percent, against the 109× already achieved, with no relationship recurring more
than seven times. A hash index per record would eat most of it back.

**The general shape of it:** git *records* history, population genetics *infers* it. So
representation and storage ideas have some chance of crossing over, and inference ideas
have none, because git has no inference to port. Having measured the representation side,
even that looks closed. The analogy transfers intuition. It does not transfer tools.

Which is a reasonable thing for an analogy to do.

---

## How anyone actually infers one

The post so far has used simulated data, where the ARG is known because it was generated
rather than inferred. Worth being explicit about what happens with real sequences.

You cannot search the space of ARGs exhaustively. It is astronomically large, and
recombination makes it worse rather than better. So the practical methods make
compromises, and the compromise each one makes is the main thing distinguishing them:

- **`tsinfer`** builds a tree sequence heuristically by matching haplotypes against each
  other. It scales to hundreds of thousands of samples and does not attempt to give you
  branch lengths in generations — dating is a separate step (`tsdate`).
- **Relate** and **ARG-Needle** infer genome-wide genealogies with time estimates, at
  biobank scale. ARG-Needle in particular was built with association testing in mind.
- **Full Bayesian samplers** exist and give you a posterior over ARGs rather than a point
  estimate, which is the statistically correct thing. They run on dozens of genomes, not
  hundreds of thousands.

That last line is the trade-off the field is currently stuck on, and the 2025 *Nature
Reviews Genetics* review says so directly: the faster methods tend to perform worse, and
only a few scale to GWAS-sized data. You can have calibration or you can have scale.

Which connects straight back to break 6. Most analyses take a point estimate from a fast
method and proceed as though the history were known. If that makes you uncomfortable,
your instinct is correct and the review shares it.

---

## Where to go from here

- **[The era of the ARG: an empiricist's guide](https://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1011110)** — the standard introduction, written for people new to genealogies.
- **Nielsen, Vaughn & Deng**, *Nature Reviews Genetics* 26:47–58 (2025) — ARG inference methods and what remains unsolved.
- **[`tskit`](https://tskit.dev) and [`msprime`](https://tskit.dev/msprime/)** — the libraries. `msprime.sim_ancestry()` gives you a real ARG in one line, and the tutorials are good.
- **[Modelling complex traits with ancestral recombination graphs](https://www.biorxiv.org/content/10.1101/2025.07.14.664631v1)** — where this connects to heritability and GWAS, and the source of the claims above about heritability not being a single number and about selection producing structure that neutral methods do not expect.

The analogy gets you to the door. Past it, the coalescent is a real probabilistic model
and worth learning properly. But knowing that you are reconstructing a deleted commit log
is a better starting point than most people get.

---

## Glossary

Terms you will hit immediately in the literature, in the order they stop being confusing.

| term | what it means |
| --- | --- |
| **locus** | a position or region in the genome. A path in the repo. |
| **allele** | one of the possible states at a locus. The content of that line. |
| **ancestral / derived** | the original state, and the state after a mutation. Before and after the diff. |
| **haplotype** | one contiguous copy of a chromosome, inherited as a unit. One checkout. |
| **SNP / variant** | a position where the sample carries more than one allele. A line that differs between checkouts. |
| **allele frequency** | the fraction of samples carrying the derived allele. The size of the clade under the branch. |
| **coalescence** | two lineages meeting in a common ancestor, viewed backward in time. |
| **TMRCA** | time to most recent common ancestor. The depth of the join. |
| **marginal / local tree** | the genealogy at one position, valid until the next breakpoint. |
| **breakpoint** | a position where recombination changes the local tree. |
| **tree sequence** | the whole set of local trees, stored as edits. What `tskit` holds. |
| **IBD** | identity by descent: a segment two samples inherited from a shared ancestor without recombination between them. |
| **linkage disequilibrium** | correlation between alleles at different loci. Tree similarity, seen through genotypes. |
| **effective population size (Ne)** | the size of an idealised population that would produce the observed amount of genetic variation. Sets the timescale of coalescence. |
| **molecular clock** | the assumption that mutations accumulate at a roughly steady rate, which is how branches get dated. |
| **infinite sites** | the assumption that each position mutates at most once in the history. Convenient, and sometimes false. |
| **homoplasy** | a shared variant that is not from a shared ancestor. The assumption above, failing. |
| **selective sweep** | a beneficial variant rising in frequency and dragging its neighbours along. |

---

## Reproducing the numbers

```bash
pip install -r measurements/requirements.txt
python measurements/headroom.py      # the compression and delta-base measurements
python measurements/linkage.py       # linkage disequilibrium vs genealogical distance
python measurements/make_example.py  # regenerates the data behind the interactive page
```

Both are single simulations with fixed seeds, stated parameters, and no tuning. They are
illustrations rather than a benchmark: the compression figures depend on sample size and
recombination rate, and a serious claim about either would need a sweep. The qualitative
results — adjacent trees are near-optimal delta bases, and `tskit`'s representation leaves
little on the table — held everywhere I looked.

## Acknowledgements

The simulations use [`msprime`](https://tskit.dev/msprime/) and
[`tskit`](https://tskit.dev), both of which are excellent and neither of which I wrote.

## License

Text CC BY 4.0, code MIT.
