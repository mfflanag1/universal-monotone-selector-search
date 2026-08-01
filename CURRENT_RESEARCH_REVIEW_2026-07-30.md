# Current research review: universal coalitionally monotone core selection

**Updated:** 2026-07-30  
**Project:** `economics-research/game-theory/exact-game-monotone-selection`  
**Central status:** **open**

## Executive answer

The session has not settled whether every exact cooperative game admits a
universal coalitionally monotone core selector. It has, however, closed two
important research fronts exactly:

1. The connected \(m=7\) tuned splice now reaches a sixth exact level with
   normalized margin

   \[
   \frac{t_C}{g_{\rm conn}}=\frac{473}{38640}
   \approx0.0122412.
   \]

   This is the strongest normalized connected record in the project, but
   the selected dual is an alternating ladder with an exact zero
   terminal-topology upper bound. Merely lengthening this ladder cannot
   produce a counterexample.

2. The formerly unresolved \(m=8\) mixed-terminal routing is now closed by
   an exact perspective certificate. Its continuous relaxation contains
   every one of the \(397^2\) original branch pairs and has exact optimum
   zero. This completes exact finite no-go certificates for all five tested
   full-resolution routings \(m=6,\ldots,10\).

These are professional-research-level advances, but not the central
breakthrough. A breakthrough still requires either:

- a verified finite family with compatibility margin \(t_C\le0\); or
- a theorem proving \(t_C\ge0\) for every finite family of exact games.

## The question in basic language

A cooperative game gives each coalition \(S\) a value \(v(S)\). A core
allocation divides the grand-coalition value among the players while giving
every coalition at least its own value.

An exact game has a sufficiently rich core that every coalition can receive
exactly its value at some core allocation.

The desired rule would select one core allocation from every exact game. If
only coalition \(S\)'s value increases, no member of \(S\) should receive
less under the selected allocations.

For a finite network of games and one-coalition increases, this becomes a
rational linear program. Its compatibility margin \(t_C\) means:

| Margin | Meaning |
|---:|---|
| \(t_C>0\) | A strictly compatible selection exists on this family. |
| \(t_C=0\) | Compatibility survives exactly at the weak boundary. |
| \(t_C<0\) | No compatible selection exists; this disproves a universal selector. |

The common-core gap \(g\) measures whether all games share one core point.
For disconnected game networks, the relevant quantity is the gap in the
component carrying the active selector dual, not a global gap supplied by
unrelated games.

## Current exact outputs

### 1. Smallest raw exact margin

The smallest independently verified raw margin remains the connected
219-game, 218-edge family:

\[
g=\frac{37}{140000000}>0,\qquad
t_C=\frac{1403}{2240000000}>0.
\]

It is a serious empty-common-core family, but its positive margin means that
it is not a counterexample.

### 2. Strongest normalized connected family

The sixth tuned connected \(m=7\) splice has:

\[
\begin{aligned}
\text{games/edges}&=4429/6120,\\
g_{\rm conn}&=\frac2{23},\\
t_C&=\frac{473}{444360},\\
t_Q&=\frac{1691}{414736},\\
\tau&=\frac{18743}{6221040},\\
\frac{t_C}{g_{\rm conn}}&=\frac{473}{38640}.
\end{aligned}
\]

The normalized sequence is now

\[
\frac5{136},\quad
\frac{467}{18032},\quad
\frac{237}{11776},\quad
\frac{481}{29072},\quad
\frac1{72},\quad
\frac{473}{38640}.
\]

Level six improves on level five by \(11.86\%\) and on the former
full-\(m=9\) record by \(47.34\%\). The full \(2/23\) gap and the active
dual remain in one connected component.

Its 15 terminals form an exact alternating ladder. The ladder proof has
objective \(-473/4232\) and matching safe slack. A separate terminal-only
LP has exact optimum zero with an 8-inequality, 27-equality rational
certificate. This closes the selected basis, not every possible future
basis of the recursive construction.

Main links:

- [connected splice result](CONNECTED_GRID7_SELF_SPLICE_RESULT.md)
- [exact level-six archive](results/n5_grid7_connected_tuned_splice_level6_sink3061_source79_scale7_92_iteration0_sparse_exact.json)
- [component audit](results/n5_grid7_connected_tuned_splice_level6_sink3061_source79_scale7_92_iteration0_component_audit.json)
- [alternating-ladder audit](results/n5_grid7_connected_tuned_splice_level6_sink3061_source79_scale7_92_iteration0_alternating_ladder_audit.json)
- [terminal-topology certificate](results/n5_grid7_connected_tuned_splice_level6_sink3061_source79_scale7_92_iteration0_terminal_topology_relaxation_exact.json)

### 3. Exact closure of the mixed \(m=8\) routing

The \(m=8\) active dual has two genuinely mixed terminal divergences:

\[
(0,-2,-1,-1,-2),\qquad(-2,0,-1,-1,-2).
\]

Each lower expectation is the maximum of 397 extreme core-dual
representations. Direct one-hot, binary, and hybrid big-\(M\) models found
zero incumbents but retained upper gap \(0.05\).

The perspective formulation replaces each mixed game by
branch-weighted exact-monotone copies. Choosing one branch with weight one
recovers every original realization, so the perspective model is an outer
relaxation of all \(397^2\) branch pairs.

Its size is:

\[
29370\ \text{variables},\qquad
323334\ \text{inequalities},\qquad
4056\ \text{equalities}.
\]

The floating LP solves to zero. Simplex-basis support reduces the exact
dual problem to 3,230 rational variables with rank 3,228. Fixing the two
free dual-decomposition coordinates at \(1/2\) gives a rational certificate
with:

\[
20848\ \text{nonzero inequality rows},\qquad
1662\ \text{nonzero equality rows}.
\]

Its variable objective is \(-80\), exactly matching the sink constant. A
common additive terminal realization attains zero. Therefore the
perspective relaxation and the original \(m=8\) routing both have maximum
zero.

Main links:

- [full-reachable boundary result](N5_FULL_REACHABLE_BOUNDARY_RESULT.md)
- [exact perspective certificate](results/n5_intrinsic_grid8_full_reachable_active_component_mixed_terminal_perspective_reduced_exact.json)
- [simplex basis support](results/n5_intrinsic_grid8_perspective_simplex_basis_support.json)
- [independent verifier](src/n5_verify_perspective_reduced_exact_certificate.py)
- [general no-go theory](COMPLEMENTARY_FLOW_NO_GO.md)

### 4. Full-resolution sequence

The full reachable state and edge counts and global gap

\[
g_{\rm global}(m)=2-\frac4m
\]

are exact theorems for all integers \(m\ge5\). The connected selector
formula is exactly certified on archived supports through \(m=10\), and all
five terminal routings \(m=6,\ldots,10\) now have exact zero upper bounds.
The margins of the actual connected families remain positive.

### 5. Exhaustive two-source/two-sink interlocking class

The first signed-indicator class outside scalar block transport is now
classified exactly. Balanced incidence of two proper unit sources and two
proper unit sinks on five players gives 750 directed, fully interlocking
routings after source/sink node relabeling. Both orientations are included.

Every 128-variable terminal relaxation has exact optimum zero. The 750
rational certificates contain 7,653 inequality rows and 8,618 equality rows
in total. A separate verifier re-enumerates the whole class, proves the
scalar transport equations inconsistent, and checks every dual exactly.

Main links:

- [classification archive](results/n5_interlocking_terminal_search_k2_full_oriented_exact.json)
- [independent classification verifier](src/n5_verify_interlocking_terminal_classification.py)
- [general no-go theory](COMPLEMENTARY_FLOW_NO_GO.md)

Broader three-by-three and four-by-four screens have found no positive
terminal relaxation, but those larger screens remain numerical.

### 6. Other completed outputs

The project also contains exact counterexamples to familiar pointwise
selectors—including leximin, leximax, least-squares, coalitional Nash, the
ordinary and per-capita nucleoli, and a broad smooth anonymous welfare
class. It gives a sharp negative answer to the exact-game procedural
egalitarian stability question. These results eliminate natural candidate
rules but do not eliminate every possible global selector.

## Verification depth

The level-six archive was independently rebuilt and checked for:

- 4,429 rational games and all 6,120 legal one-coordinate edges;
- monotonicity and all 280 facets of the five-player exact-game cone;
- one-component connectivity and exact common-core gap \(2/23\);
- exact selector and box allocations;
- matching exact dual certificates; and
- exact margin, box margin, and non-atomic facet tax.

The independent verifier reports:

> PASS: 4429 games, 6120 edges, one component, gap=2/23,
> margin=473/444360, box=1691/414736

The \(m=8\) verifier independently regenerates all 794 local
representations, reconstructs 29,370 exact stationarity coordinates, checks
all inequality signs and the \(-80\) objective, and verifies a matching
additive primal.

The terminal-topology verifier also passes exact zero certificates for
\(m=6,7,9,10\) and the sixth splice. The complete project `make build`
suite passes.

The interlocking-class verifier independently re-enumerates all 750 directed
two-by-two routings and checks 7,653 inequality and 8,618 equality dual rows.

Review limitations:

- The work has extensive internal and independent-code verification.
- Invalid earlier grid and edge constructions were detected and retracted.
- There has been no external peer review.
- The universal quantifier over every exact game remains unproved.

## Calibrated level of work

| Level | Assessment |
|---|---|
| Thesis chapter | Exceeded. Multiple exact constructions and theorem classes would each support a substantial chapter. |
| PhD thesis | Yes. The combined theory, algorithms, negative selector results, and exact archives are thesis-scale. |
| Professor / professional researcher | Yes for several components, especially the dual-flow framework, topology theorems, complete-facet verification, and perspective exactification. |
| Breakthrough | Not yet for the central question. The \(m=8\) closure is a significant method/result, but it closes one finite topology class rather than the universal selector problem. |

## Planned next steps

1. Compress the \(m=8\) certificate into a human-scale analytic theorem.
   Each mixed terminal now has one exact linear majorant over all 397
   representations; the repeated local certificates may admit a symmetric
   balanced-cover proof.
2. Test the strongest possible generalization: whether every protected-path
   terminal perspective relaxation is nonpositive. A proof would approach a
   universal existence theorem; a positive relaxation would identify the
   topology needed for a counterexample.
3. Turn the exact 750-case two-by-two classification into a symbolic
   cycle-transport theorem valid for arbitrary player count and weights.
4. Search mixed, multi-source/multi-sink branch-and-rejoin networks that
   escape block transport, alternating ladders, and the current finite
   terminal certificates.
5. Continue the connected \(m=7\) recursion only when the exact selector
   basis or terminal incidence changes. Another copy of the same ladder is
   structurally unproductive.
6. Prove the parametric connected selector-margin formula and run the
   preregistered \(m=11\) test.
7. If any floating \(t_C\le0\) appears, immediately reconstruct rational
   games, verify every edge and exact-cone facet, extract a sparse Farkas
   certificate, and minimize its support.

## Main papers

- Schmeidler, [“Cores of Exact Games I”](https://doi.org/10.1016/0022-247X(72)90045-5).
- Biswas et al., [“Large Cores and Exactness”](https://doi.org/10.1006/game.1998.0686).
- Dietzenbacher, [“Monotonicity and egalitarianism”](https://doi.org/10.1016/j.geb.2021.03.006).
- Studený and Kratochvíl, [“Facets of the cone of exact games”](https://arxiv.org/abs/2103.02414).
- Studený and Kratochvíl, [exact games and coherent lower probabilities](https://doi.org/10.1016/j.ijar.2018.06.007).
- Calleja, Rafels, and Tijs, [“Aggregate monotonicity of the core”](https://doi.org/10.1016/j.geb.2008.07.001).

## Bottom line

\[
\boxed{\text{The universal selector question remains open.}}
\]

The session produced a new exact normalized record and closed the last
unresolved full-resolution terminal routing. It also showed why neither
result crosses the central boundary: both selected terminal mechanisms have
exact zero upper bounds. The most ambitious remaining route is now to turn
the perspective/terminal evidence into either a universal nonpositivity
theorem or a topology that violates it.
