# Sharp mixed-divergence premiums on five-player exact games

**Date:** 2026-08-01  
**Central selector question:** open

## Result

Let \(v\) be a normalized monotone exact game on five players and let

\[
\underline E_v(d)=\min_{x\in C(v)}d\cdot x.
\]

Suppose \(d\in\{-1,0,1\}^5\) uses all three levels.  Write

\[
L=\{i:d_i=-1\},\qquad M=\{i:d_i=0\},\qquad H=\{i:d_i=1\}.
\]

The nested-level-set or Choquet lower bound is

\[
\mathcal C_v(d)=-v(N)+v(N\setminus L)+v(H).
\]

The sharp maximum of the non-additivity premium

\[
\Delta_v(d)=\underline E_v(d)-\mathcal C_v(d)
\]

over the entire normalized monotone exact-game cone depends only on the
level multiplicities:

| \((|L|,|M|,|H|)\) | Sharp \(\max_v\Delta_v(d)\) | Extreme branches |
|---:|---:|---:|
| \((1,1,3)\) | \(1/3\) | 409 |
| \((1,2,2)\) | \(1/2\) | 397 |
| \((1,3,1)\) | \(1/2\) | 310 |
| \((2,1,2)\) | \(0\) | 513 |
| \((2,2,1)\) | \(1/2\) | 397 |
| \((3,1,1)\) | \(1/3\) | 409 |

Thus the Choquet replacement is false in five of the six three-level types.
The failure can be as large as one half of the grand worth.  This rules out
the simplest proposed universal proof, which would replace each mixed
divergence by its nested indicator decomposition.

The middle split is exceptional and exact:

> If two players have the low value, one player has the middle value, and
> two players have the high value, then
> \(\underline E_v(d)=\mathcal C_v(d)\) for every five-player exact game.

No game-monotonicity inequalities are needed for this zero-premium identity.
Equivalently, if \(A\subset U\subset N\), \(|A|=2\), and \(|U|=3\), then
every five-player exact game has a core allocation \(x\) satisfying

\[
x(A)=v(A),\qquad x(U)=v(U).
\]

Indeed, minimizing \(\mathbf1_A+\mathbf1_U\) has lower bound
\(v(A)+v(U)\), and the 2-1-2 identity shows this bound is attained.  The
same common tight point proves the identity for arbitrary three levels
\(a<b<c\), not only equally spaced levels.

The next degree class is now classified as well.  If \(d\) has four equally
spaced levels \(0,1,2,3\), all four possible positive multiplicity patterns
have the same sharp premium:

| level multiplicities | Sharp \(\max_v\Delta_v(d)\) | Extreme branches |
|---:|---:|---:|
| \((1,1,1,2)\) | \(1/2\) | 599 |
| \((1,1,2,1)\) | \(1/2\) | 549 |
| \((1,2,1,1)\) | \(1/2\) | 549 |
| \((2,1,1,1)\) | \(1/2\) | 599 |

These are exactly the new local types created by a unit-flow terminal of
degree three on five players.  The complete four-level certificate covers
2,296 additional lower-expectation branches.

## Exact proof

For fixed \(d\), the core LP dual writes every extreme lower-expectation
branch as

\[
\beta v(N)+\sum_S\alpha_Sv(S),\qquad
d=\beta\mathbf1_N+\sum_S\alpha_S\mathbf1_S,quad \alpha_S\ge0.
\]

There are at most four positive proper-coalition coefficients in an extreme
branch.  The verifier enumerates every nonsingular choice and checks the
representation over the rationals.

For every branch it then gives an exact conic identity proving

\[
\beta v(N)+\sum_S\alpha_Sv(S)
\le \mathcal C_v(d)+M v(N),
\]

where \(M\) is the value in the table.  Each right-minus-left coefficient
vector is reconstructed as a nonnegative rational combination of the 280
five-player exact-cone facets, the elementary game-monotonicity inequalities,
and the empty-set equality.  Matching rational witness games and branches
attain every stated bound.

The complete certificate covers 2,435 branches.  The special 2-1-2 archive
covers its 513 branches with only exact-cone facets: it uses neither game
monotonicity nor normalization beyond the usual additive lineality.

Main files:

- `results/n5_three_level_premium_exact.json`
- `n5_three_level_premium_certificate.py`
- `results/n5_middle_split_identity_exact.json`
- `n5_middle_split_identity_certificate.py`
- `n5_mixed_divergence_premium.py`
- `results/n5_four_level_premium_exact.json`
- `n5_four_level_premium.py`

The exact verification commands are:

```bash
python n5_three_level_premium_certificate.py \
  results/n5_three_level_premium_exact.json
python n5_middle_split_identity_certificate.py \
  results/n5_middle_split_identity_exact.json
python n5_three_level_premium_certificate.py \
  results/n5_four_level_premium_exact.json
```

## Accessibility of the sharp witnesses

The sharp witness games are valid monotone exact games, but most lie on very
thin exact-cone faces.  The exact coordinate-move audit finds:

- the three \(1/2\)-premium witnesses have no nonzero legal one-coordinate
  move in either direction;
- the \((1,1,3)\) witness has upward moves for only one of its three positive
  divergence players;
- the \((3,1,1)\) witness has no upward move containing its positive player.

This does **not** make the premiums irrelevant: convexly mixing a witness
with a sufficiently interior exact game preserves nearly all of the premium
and gives small legal moves.  It does show why optimizing the local premium
alone lands on unusable boundary games and why an actual counterexample must
optimize premium and protected-path geometry jointly.

The audit is independently reproducible with:

```bash
python n5_premium_witness_bump_audit.py \
  results/n5_three_level_premium_exact.json \
  --output results/n5_mixed_premium_witness_bump_audit.json
```

## New topology no-go

The sharp positive premiums might suggest putting mixed games between
indicator sources and sinks.  A new protected block-flow forest theorem in
`COMPLEMENTARY_FLOW_NO_GO.md` shows that no acyclic version can work.  After
grouping equal-mass protected player paths into block arcs, any network whose
underlying undirected multigraph is a forest has nonpositive total lower
expectation, even with arbitrarily many mixed transshipment nodes.

The proof repeatedly removes a leaf.  Its signed-indicator lower expectation
can be moved to its neighbor using protected-path invariance and directed
worth monotonicity.  Superadditivity then absorbs that indicator into the
neighbor's mixed divergence.  The tree eventually contracts to one node with
zero divergence.

The single mixed-hub theorem is the star case.  In basic terms, one core
allocation at the hub simultaneously upper-bounds all source rewards and
lower-bounds all sink costs.  Choosing the hub allocation that minimizes its
mixed divergence makes those terms cancel.

Therefore a counterexample now has a sharper necessary structure:

1. every protected block-path decomposition must retain an undirected cycle,
   counting distinct parallel blocks as distinct multigraph edges;
2. at least one mixed node must avoid the zero-premium 2-1-2 type; and
3. the positive local premiums must survive the exact-cone and
   protected-path relations between those nodes.

This precisely matches the previously difficult \(m=8\) topology, which had
two mixed terminals of the \((2,2,1)\) type.  Its exact perspective
certificate still gives global upper bound zero, so a positive local premium
is necessary but not sufficient.

The first cycle with four simultaneously maximal local premiums is now also
closed exactly.  The Johnson-square blocks

\[
012,\quad013,\quad024,\quad034
\]

produce two \((1,2,2)\) sources and two \((2,2,1)\) sinks, each with premium
\(1/2\).  Protected equalities collapse the four-game problem to a two-game
cross-cosingleton swap.  An exhaustive exact audit certifies all
\(397^2=157{,}609\) affine branch pairs with upper bound zero.  See
`N5_JOHNSON_SQUARE_NO_GO_RESULT.md` and
`results/n5_cross_cosingleton_branch_sweep_exact.json`.

More generally, every alternating protected cycle is harmless in arbitrary
dimension, for arbitrary overlaps, arbitrary positive edge weights, and
arbitrary grand worths.  Pairing each source with one sink and applying
exactness at the unmatched block and the complement of the matched block
makes the worth terms telescope around the cycle; grand-worth differences
are individually nonpositive.  Together with leaf contraction from the
forest theorem, this closes every protected block-flow pseudoforest.  See
`ALTERNATING_CYCLE_NO_GO_RESULT.md`.

The first irreducible cycle-rank-two class is now closed at common grand
worth as well. There are exactly 20 five-player \(K_{2,3}\) symmetry orbits
in which both degree-three sources have four levels, every pair within a
source row has premium \(1/2\), and all three sinks have premium \(1/2\).
An exact balanced-envelope classification proves that all five games share a
core point in every orbit. The audit covers 66,108 rational LP duals with zero
failures. Thus even simultaneous maximal local premiums at all five terminals
do not survive the protected relations. See `N5_K23_COMMON_CORE_RESULT.md`.

## Literature position

The calculation uses the complete five-player facet description associated
with Studený and Kratochvíl, *Facets of the cone of exact games*:
<https://arxiv.org/abs/2103.02414>.  That paper supplies the general
semi-balanced facet theory and the low-dimensional facet catalogue.  The
specific sharp premium table and the 2-1-2 common-tight-face consequence
were not located in the targeted search and should be treated as new pending
external review.

A June 2026 paper by Hokari and Ishida gives a closely related but different
five-player incompatibility for population-monotonic allocation schemes:
<https://doi.org/10.1007/s00182-026-00990-6>.  Its changes between subgames
alter many worth coordinates at once and do not satisfy protected-path
invariance, so it does not directly settle the present one-coordinate-bump
problem.

The paper's appendix values make the minimal obstruction transparent.  Let
players be indexed \(1,\ldots,5\).  In the two-player subgame on \(\{3,4\}\),
efficiency requires

\[
x^{34}_3+x^{34}_4=6.
\]

In the four-player subgame on \(1234\), the coalition-\(123\) core constraint
and grand worth 11 force \(x^{1234}_4\le2\).  Population monotonicity along
\(34\subset234\subset1234\) therefore gives \(x^{34}_4\le2\).  In the
three-player subgame on \(345\), the coalition-\(45\) constraint and grand
worth 9 force \(x^{345}_3\le3\); monotonicity along \(34\subset345\) gives
\(x^{34}_3\le3\).  Hence the two-player total is at most 5, contradicting 6.

This is exactly the complementary-fork shape that repeatedly appeared in
the present project.  It cannot survive a fixed-player exact lift: two
protected paths from one source to the cap games preserve the complementary
coalition worths at the source, where balancedness makes the two cap bounds
sum to at least the source budget.  The external counterexample is therefore
useful as a negative control and confirms, rather than bypasses, the
protected-fork no-go theorem.

## Consequence for the central question

The universal coalitionally monotone selector is neither proved nor
disproved here.  The result closes the local analytic subproblem exactly and
eliminates two overly broad proof/counterexample routes:

- arbitrary mixed lower expectations cannot all be replaced by Choquet
  values; and
- one mixed hub can never exploit the positive premium.

The remaining target is an unequal-grand or non-sharp block-flow component of
cycle rank at least two, a more general rank-two orientation, or a six-player
topology outside the five-player premium classification. Unequal grand worths
do not rescue pseudoforests.
The next proof step is to determine whether the exact facet corrections at
degree-three terminals can always be charged to protected path slack.
Failure of that charging identity would provide the algebraic template for a
counterexample.
