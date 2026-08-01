# Two-chart rolling-triple theorem

**Date:** 2026-07-30  
**Central selector question:** open

## Theorem

For every real parameter \(L\in[5,7]\), there is a finite family of 43
five-player exact games joined by 136 legal one-coalition increases such
that

\[
v_L(N)=6+\frac4L,\qquad
\operatorname{gap}(L)=2-\frac4L,\qquad
t_C(L)=\frac{2}{33L}.
\]

Here \(\operatorname{gap}(L)\) is the exact minimum budget required by a
common allocation minus \(v_L(N)\), and \(t_C(L)\) is the sharp max-min
coalitional-monotonicity margin.

In particular, every family has an empty common core and a strictly positive
selector margin. This is an exact continuum of increasingly difficult
compatible families, not a counterexample to a universal selector.

The empty common core in this statement is global across all archived games.
A later weak-component audit finds that the selector components at
\(L=5,6,7\) have gaps \(-4/5,-2/3,-4/7\), respectively. The envelope games
that create the positive global gap are disconnected from the active flow.
Thus the theorem is an exact parametric-LP result, but not a continuum of
connected empty-core obstructions.

## Construction

Three independently verified integer-resolution archives are used:

| Rolling resolutions | Common grand worth | Games / edges | Gap | \(t_C\) |
|---|---:|---:|---:|---:|
| \(5,6,7\) | \(34/5\) | 43 / 136 | \(6/5\) | \(2/165\) |
| \(6,7,8\) | \(20/3\) | 43 / 136 | \(4/3\) | \(1/99\) |
| \(7,8,9\) | \(46/7\) | 43 / 136 | \(10/7\) | \(2/231\) |

The first pair is isomorphic as a directed graph with every edge labeled by
its changed coalition. Affine interpolation of games and allocations in
\(1/L\) gives a chart on \(5\le L\le6\). The second pair is likewise
label-isomorphic and gives a chart on \(6\le L\le7\). The node mappings
change at \(L=6\), where the active exact-cone regime changes.

The initial \(m=7,8,9\) dual support had 44 games and 138 edges. A separate
sparse-dual solve proved that one game and two edges were redundant, yielding
the required 43/136 graph without changing either exact optimum.

## Exact verification

For each chart, `n5_parametric_rolling_triple_verify.py` checks symbolically
over the entire real interval:

1. game monotonicity;
2. all 280 facets of the complete five-player exact-game cone;
3. strict positivity and one-coordinate legality of all 136 edges;
4. a rational parametric core allocation attaining \(2/(33L)\);
5. a fixed rational dual with the matching objective;
6. a common allocation of budget \(8\); and
7. a balanced lower certificate proving that budget \(8\) is minimal.

The formulas then follow exactly:

\[
8-\left(6+\frac4L\right)=2-\frac4L
\]

for the common-core gap, while the matching primal and dual give
\(t_C(L)=2/(33L)\).

An additional exhaustive domain calculation checks 18,736 affine
inequalities per chart. It proves that the archived first chart is feasible
exactly on \([5,6]\) and the archived second chart exactly on \([6,7]\).
Core-allocation constraints bind at both outer endpoints; exact-cone facets
also block the first chart at \(L=5\). Thus neither chart can simply be
extrapolated. Extending the theorem requires a new \(m=8,9,10\) chart.

## Interpretation

This theorem closes two possible loopholes in the discrete evidence:

- the values at \(L=5,6,7\) are not unrelated rational coincidences; and
- positivity is not caused by integrality of the resolution parameter.

It does not prove the apparent pattern for all \(L\), and it does not settle
the universal-selector question. The margin remains positive throughout the
proved domain and tends toward zero only if further charts continue the law.
It must not be used as a normalized counterexample frontier because its
active selector components admit common-core selections.

## Reproducible outputs

- `/private/tmp/n5_parametric_rolling_triple_verified.json`
- `/private/tmp/n5_parametric_rolling_triple_6_7_verified.json`
- `/private/tmp/n5_parametric_rolling_triple_5_6_domain.json`
- `/private/tmp/n5_parametric_rolling_triple_6_7_domain.json`
- `/private/tmp/n5_intrinsic_grid789_common_grand46_7_permutation_pair_compact_exact.json`

The `/private/tmp` paths are working outputs. Final copies belong under the
economics research project’s `results/` directory.
