# Exact common-core classification of sharp five-player K2,3 topologies

**Status:** exact finite theorem  
**Scope:** all irreducible maximally sharp five-player \(K_{2,3}\) block flows
with common grand worth  
**Central selector question:** still open

## Result

Let two source games feed three sink games through six nonempty proper
protected blocks. Consider the maximally sharp irreducible class defined by:

- each degree-three source divergence has all four incidence levels;
- every pair of blocks within either source row has sharp local premium
  \(1/2\);
- every degree-two sink divergence has sharp local premium \(1/2\); and
- the union of the six blocks is \(N\), while their intersection is empty.

After normalizing the first source row by sink relabeling, there are exactly
4,440 configurations. Quotienting by player permutations, source exchange,
and sink permutations leaves exactly 20 symmetry orbits.

For every one of those 20 orbits, assume the five games are normalized
monotone exact games with a common grand worth and impose every worth equality
and inequality forced by the protected paths. Then

\[
\boxed{\bigcap_{k=0}^4 C(v^k)\ne\varnothing.}
\]

Consequently every terminal-flow lower-expectation objective is
nonpositive. If the terminal divergences are \(d_k\) and sum to zero, a
common core point \(x\) gives

\[
\sum_k\underline E_{v^k}(d_k)
\le \sum_k d_k\cdot x=0.
\]

This is stronger than closing the unit-flow divergence used to discover each
topology: it closes every zero-sum divergence supported on the same five games
and protected relations.

One maximally unconstrained representative, used as the initial target, has
the six protected blocks

\[
\begin{pmatrix}
\{2\}&\{0,2,4\}&\{0,1,2\}\\
\{0,2,3\}&\{0\}&\{0,1,3\}
\end{pmatrix},
\]

or masks

\[
\begin{pmatrix}4&21&7\\13&1&11\end{pmatrix}.
\]

## Why this was the first genuinely unresolved target

The protected block-flow forest and alternating-cycle theorems close every
pseudoforest. The displayed \(K_{2,3}\) has six edges and five vertices, so
its cycle rank is two. Each degree-three source has all four possible
incidence levels, every pair within either source row has sharp five-player
premium \(1/2\), and every degree-two sink also has sharp premium \(1/2\).
In particular, the zero-premium 2-1-2 common-tight identity does not remove
one of its two cycles.

## Exact finite proof

Because each sink game dominates both source games coordinatewise and the
grand worth is common, each sink core is contained in both source cores.
It is therefore enough to prove that the three sink cores intersect.

Define their pointwise worth envelope

\[
w(S)=\max\{v^{q_0}(S),v^{q_1}(S),v^{q_2}(S)\}.
\]

Then

\[
C(w)=C(v^{q_0})\cap C(v^{q_1})\cap C(v^{q_2}).
\]

For the displayed hardest representative, protected-path invariance forces
the three sink worths to agree outside only six proper coalitions, with masks

\[
7,\ 13,\ 15,\ 21,\ 23,\ 29.
\]

Balancedness of \(w\) can be checked on the extreme points of the balanced-
weight polytope

\[
\lambda\ge0,
\qquad
\sum_S\lambda_S\mathbf1_S=\mathbf1_N.
\]

The verifier enumerates these vertices using integer determinants and
Cramer's rule. There are exactly 1,291. Every vertex has support at most
five. On coordinates where the sinks agree, the maximum in \(w\) is fixed;
on a remaining support coordinate, it is one of three sink worths. Depending
on the orbit, only two through six sink coordinates can differ.

For every case, the verifier maximizes the selected balanced worth over the
full five-game protected topology, including all 280 exact-game facets per
node and all elementary game-monotonicity inequalities. Every optimum is at
most the grand worth. Across all 20 orbits there are 66,108 envelope cases.
All 66,108 LP duals reconstruct exactly over the rationals, with zero
failures. A common additive realization attains equality, so every exact
optimum is one and every common-core gap is exactly zero.

The independent coverage verifier re-enumerates all 4,440 sink-normalized
configurations, reconstructs the 20 symmetry orbits, canonicalizes every
certificate, rejects duplicates or omissions, and checks the exact status and
case count of every orbit.

## Reproduction

```bash
python n5_k23_common_core_sweep.py \
  --output results/n5_k23_all_pairs_sharp_common_core_exact.json
python n5_verify_k23_sharp_common_core_classification.py \
  --output results/n5_k23_all_sharp_common_core_classification_exact.json
```

Main files:

- `n5_k23_common_core_sweep.py`
- `n5_verify_k23_sharp_common_core_classification.py`
- `results/n5_k23_all_pairs_sharp_common_core_exact.json`
- `results/n5_k23_all_sharp_common_core_classification_exact.json`
- `results/n5_balanced_vertices_exact.json`
- `n5_k23_all_pairs_sharp_terminal_report.json`

## Boundary of the theorem

This result is exact but deliberately bounded. It does not yet cover unequal
grand worths, non-sharp five-player \(K_{2,3}\) configurations, arbitrary
cycle-rank-two block flows, or the universal selector question. Numerical
sweeps of the same 20 sharp orbits with unequal grand worths found no positive
obstruction; that unequal-grand sweep is evidence, not part of this theorem.

The common-core method genuinely stops outside the classified sharp class.
For masks

\[
\begin{pmatrix}4&15&8\\4&23&8\end{pmatrix},
\]

the exact unanimity-game realization with masks

\[
12,\ 12,\ 4,\ 12,\ 8
\]

has common-core budget gap exactly one. The singleton-unanimity sinks force
both player 2 and player 3 to receive one against grand worth one. Yet the
terminal lower-expectation objective at this realization is \(-3\), not
positive, and a deep direct topology search still maximizes at zero. Thus an
empty common core is not enough; the empty-core separation must align with the
protected terminal divergence. See
`results/n5_k23_empty_common_core_witness_exact.json` and
`n5_verify_k23_empty_common_core_witness.py`.
