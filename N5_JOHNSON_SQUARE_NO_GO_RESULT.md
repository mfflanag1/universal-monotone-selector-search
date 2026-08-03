# Exact five-player Johnson-square no-go theorem

Date: 2026-08-01

## Result

The first four-terminal protected-flow cycle in which all four local
divergences attain a sharp positive five-player nonadditivity-premium type is
nevertheless harmless.  Its total lower-expectation objective is always
nonpositive on normalized monotone exact five-player games.

This closes a genuine cycle beyond the protected block-flow forest theorem.
It does **not** settle the universal selector question.

Write players as (0,1,2,3,4).  The four protected blocks are

\[
B_{00}=012,\qquad B_{01}=013,\qquad
B_{10}=024,\qquad B_{11}=034,
\]

on the directed (K_{2,2}) topology

\[
s_0\longrightarrow t_2,t_3,
\qquad
s_1\longrightarrow t_2,t_3.
\]

The terminal divergences are

\[
\begin{aligned}
d_{s_0}&=(2,2,1,1,0),&
d_{s_1}&=(2,0,1,1,2),\\
d_{t_2}&=(-2,-1,-2,0,-1),&
d_{t_3}&=(-2,-1,0,-2,-1).
\end{aligned}
\]

The two source types and two sink types each have sharp local premium
(1/2).  Local premium bounds alone therefore leave a nominal total premium
of (2); the theorem is a genuinely coupled result.

## Cross-cosingleton swap lemma

Let (p,q) be normalized monotone exact five-player games that agree on
every coalition except possibly

\[
A=N\setminus\{2\},\qquad B=N\setminus\{3\},
\]

and satisfy

\[
p(A)\le q(A),\qquad q(B)\le p(B).
\]

For

\[
e=(2,1,0,2,1),
\]

the following exact inequality holds:

\[
\boxed{\underline E_p(e)+\underline E_q(-e)\le0.}
\]

The lemma says that the ranges of the linear statistic (e\cdot x) over the
two cores cannot cross in the wrong order:

\[
\min_{x\in C(p)}e\cdot x
\le
\max_{y\in C(q)}e\cdot y.
\]

### Exact finite certificate

For either divergence, the core LP dual has 397 extreme representations.
Every representation consists of the free grand-coalition constant and at
most four positive proper-coalition indicators.  The enumerator checks its
defining equations exactly after discovering each nonsingular (0/1) basis.

For every one of the

\[
397^2=157{,}609
\]

pairs of extreme representations, the certificate program maximizes the
corresponding affine sum over the complete paired five-player polytope:

- all 280 exact-cone facets at each node;
- all game-monotonicity inequalities;
- equality away from (A,B); and
- the two directed cross inequalities at (A,B).

Every affine maximum has an exact rational dual upper bound of zero.  The
program reconstructs the sparse dual multipliers as `Fraction` objects and
checks, without tolerance:

1. nonpositivity of every inequality multiplier in the minimization dual;
2. exact stationarity in all 64 game-worth coordinates; and
3. exact dual objective zero.

The floating discovery pass had largest residual
(2.8892443992845074\times10^{-12}); this number is not used as proof.  The
rational dual identities are the proof.

Run

```sh
python n5_cross_cosingleton_branch_sweep.py \
  --exact-audit \
  --output results/n5_cross_cosingleton_branch_sweep_exact.json
```

The run succeeds only if all 157,609 exact dual reconstructions pass.

## Reduction of the Johnson square to the swap lemma

For a fixed coalition (S), an edge of the block square forces equality of
its endpoint worths unless its protected block is contained in (S).  The
inactive edges connect all four terminals except in four proper-coalition
cases.  Consequently the four games agree everywhere except:

| coalition | only possible exceptional terminal | direction |
|---|---:|---:|
| (N\setminus\{4\}=0123) | (s_0) | lower |
| (N\setminus\{1\}=0234) | (s_1) | lower |
| (N\setminus\{3\}=0124) | (t_2) | higher |
| (N\setminus\{2\}=0134) | (t_3) | higher |

In particular, the two sink games (t_2,t_3) satisfy the
cross-cosingleton lemma: (t_2) may be higher only at
(N\setminus\{3\}), and (t_3) may be higher only at
(N\setminus\{2\}).

Lower expectation is monotone in the game worths and superadditive in its
divergence argument.  Since (s_0\le t_2), (s_1\le t_2), and

\[
d_{s_0}+d_{s_1}+d_{t_2}=-d_{t_3}=e,
\]

we obtain

\[
\begin{aligned}
\sum_k\underline E_k(d_k)
&\le
\underline E_{t_2}(d_{s_0})
+\underline E_{t_2}(d_{s_1})
+\underline E_{t_2}(d_{t_2})
+\underline E_{t_3}(d_{t_3})\\
&\le
\underline E_{t_2}(e)+\underline E_{t_3}(-e)\\
&\le0.
\end{aligned}
\]

This proves the Johnson-square no-go theorem.

## No-common-player square audit

Removing the common player leaves 68 permutation/row/column orbits of
five-player unit-flow block squares in which all four terminal types have the
maximal local premium \(1/2\).  A continuous coordinate-ascent audit used
eight random/structured starts and twelve branch updates on every orbit.  No
positive objective was found; every best value was zero to numerical
tolerance.

A stricter 22-orbit subfamily requires four distinct incomparable size-two or
size-three blocks, full five-player support, and no common player.  A deeper
run with 24 random/structured starts and 20 updates per orbit again found no
positive objective.

These two audits are search evidence, not exact no-go certificates.  Their
archives are

- `results/n5_all_max_premium_block_square_batch.json`; and
- `results/n5_dangerous_block_square_batch_deep.json`.

They indicate that simply deleting the Johnson square's common player is not
enough.  The next topology should also escape the complete-bipartite
two-sink interval reduction.

## Research consequence

The first cyclic topology that simultaneously realizes four maximally
dangerous local five-player types still cannot obstruct compatibility.  The
failure is not explained by a weak local terminal: it comes from exact
cross-game coupling, which forces the four-game problem through a two-game
cross-cosingleton interval-overlap inequality.

The next search should therefore remove the shared-player geometry, use a
cycle whose inactive-edge graph leaves more than one private worth coordinate
per terminal, or move to six players where the five-player facet and premium
classification no longer applies.
