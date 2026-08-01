# Full reachable-component boundary in the five-player exact cone

**Date:** 2026-07-30  
**Central selector question:** open

## Result

The intrinsically five-player exactification grid has now been solved using
every exact state reachable below its sharp grand cap, rather than one
chosen coordinate path. Five independently verified exact resolutions give:

| \(m\) | Full games / edges | Active component | Component gap | \(t_C\) | \(t_Q\) |
|---:|---:|---:|---:|---:|---:|
| 5 | 13,755 / 64,916 | 195 / 291 | \(-2/15\) | \(17/620\) | \(7/135\) |
| 6 | 18,963 / 90,448 | 394 / 584 | \(0\) | \(13/564\) | \(2/45\) |
| 7 | 24,171 / 115,980 | 370 / 510 | \(2/7\) | \(5/252\) | \(3/77\) |
| 8 | 29,379 / 141,512 | 494 / 679 | \(1/2\) | \(11/632\) | \(3/76\) |
| 9 | 34,587 / 167,044 | 897 / 1,409 | \(2/3\) | \(53/3420\) | \(13/396\) |
| 10 | 39,795 / 192,576 | 1,275 / 2,074 | \(4/5\) | \(31/2220\) | \(3/110\) |

The core margins satisfy

\[
\boxed{
t_C(m)=\frac{9m-28}{4m(16m-49)}
}
\qquad(m=5,6,7,8,9,10).
\]

The \(m=8\), \(m=9\), and \(m=10\) values were predicted before their solves
and all passed out of sample. The \(m=10\) result is independently exactified
on a connected 1,275-game support and verified with

\[
t_C=\frac{31}{2220},\qquad
t_Q=\frac3{110},\qquad
g_{\rm conn}=\frac45.
\]
The full graph counts also satisfy

\[
\#V_m=5208m-12285,\qquad
\#E_m=25532m-62744.
\]

### Exact state-level recurrence audit

The count law is not merely a regression on five totals.
n5_full_m_recurrence_audit.py reconstructs the ten-dimensional integer
state of every game on both branches. For every archived
\(m=5,\ldots,10\), each branch's complete reachable set is exactly

\[
P\;\cup\!
\bigcup_{k=0}^{m-6}(R+k\mathbf1)
\;\cup\;(T+(m-5)\mathbf1),
\]

where the three fixed rationally reconstructed templates have

\[
|P|=3688,\qquad |R|=2604,\qquad |T|=3190.
\]

The repeated template occupies ten consecutive depth layers with histogram

\[
(305,304,299,265,247,226,224,220,241,273),
\]

whose sum is 2,604. The two branches share only the origin, so this exact
finite audit implies the displayed affine vertex law on every archived
resolution. The induced-edge counts simultaneously pass the affine law with
increment 25,532.

This finite set identity supplies the induction templates for the all-\(m\)
theorem below. Its three symbolic obligations are exact-cone feasibility,
reachability, and exclusion of every outward boundary successor on
\(P,R,T\).

All three combinatorial proof obligations are now complete.
n5_full_m_template_verify.py checks 6,827,040 parametric inequalities. For
every real \(m\ge5\), it proves monotonicity and all 280 exact-cone facets
for every game in

\[
P\cup\bigcup_{k=0}^{m-6}(R+k\mathbf1)
\cup(T+(m-5)\mathbf1).
\]

For every integer \(m\ge5\), it also proves inductively that the unit-bump
graph reaches every displayed state from the origin. Counting all legal
one-coordinate increases—not only adjacent grid steps—gives 32,458 edges
per branch at \(m=5\) and exactly 12,766 more per added \(R\) block.
An exact interaction audit rules out every nonlocal template edge: apart
from within-template edges, the only cross-block interactions are
\(P\) to the first \(R\), adjacent \(R\) blocks, and the last \(R\) to
the shifted \(T\). The two branches change disjoint proper-coalition sets
in strictly positive directions, so two non-origin branch states can never
differ in only one coalition; there are no hidden cross-branch edges.
Thus the constant edge increment is symbolic, not an extrapolation from the
first two resolutions.
Consequently this explicit all-\(m\) exact connected family has

\[
\boxed{
\#V_m=5208m-12285,\qquad
\#E_m=25532m-62744
}
\quad(m\ge5).
\]

The same verifier proves the global common-core law for this all-\(m\)
family. The fixed budget-\(8\) allocation

\[
(1,1,3,3,0)
\]

dominates every proper-coalition worth on all three parametric templates.
Conversely, explicit terminal-template games attain the five singleton
values \(1,1,3,3,0\). The singleton partition is therefore a balanced
lower certificate of value \(8\). Since \(v(N)=6+4/m\),

\[
\boxed{g_{\rm global}(m)=2-\frac4m}
\qquad(m\ge5).
\]

Finally, n5_full_m_exclusion_verify.py checks every outward
unit-coordinate boundary move from \(P\), an interior or final copy of
\(R\), and the shifted \(T\). Every candidate outside the displayed
templates has a monotonicity or exact-cone row whose numerator is strictly
negative over its entire parameter domain. The exact \(m=5,6\) enumerations
are the base cases. Induction therefore proves the set identity, and hence
the vertex and edge laws, for the *full* reachable component at every
integer \(m\ge5\); there are no additional reachable states outside the
templates.

What remains open in this family is the selector theorem: a parametric
primal-dual proof of the margin formula.

A block-classification audit of the exact \(m=6,7,8,9\) active components
rules out the easiest such proof. Their selected \(P\), \(R_k\), and \(T\)
sets change substantially with \(m\). Instantiating the union of every base
state used at those four resolutions in all five \(m=10\) repeated blocks
produces a connected 3,182-game, 7,406-edge subfamily with gap \(4/5\), but
its numerical margin is \(0.016744186\ldots>31/2220\). The predicted \(m=10\) optimum
therefore requires new active states rather than a stationary sparse basis.

The exact core dual has normalization \(8(16m-49)\) and \(12m-40\) core
rows. At \(m=9\), its selected sparse basis has 123 nonzero-divergence nodes
rather than the \(22m-79=119\) continuation of the \(m\le8\) bases. The
margin formula therefore survives a dual-topology breakpoint. It has the
positive asymptotic law

\[
t_C(m)\sim\frac9{64m}.
\]

At \(m=10\), the exact dual scale is 888 and all 152 nonzero-divergence
nodes are signed indicators. There are 81 sources and 71 sinks. Exact
coalition-type transport LPs are infeasible in both the fork and join
orientations, so this is a genuinely interlocking indicator network outside
the block-transport no-go theorem. Its positive margin shows that escaping
that theorem is necessary but not sufficient for a counterexample.

The terminal graph is nevertheless much smaller topologically than the raw
terminal count suggests. Join a source role to a sink role when at least one
protected player path connects them, then repeatedly delete degree-one
vertices. The exact kernel audit gives:

| \(m\) | Terminal nodes | Source–sink pairs | Cycle rank | Two-core \(V/E\) | Branches |
|---:|---:|---:|---:|---:|---:|
| 6 | 54 | 55 | 2 | 25/26 | 1 |
| 7 | 73 | 74 | 2 | 33/34 | 1 |
| 8 | 99 | 102 | 4 | 45/48 | 3 |
| 9 | 152 | 163 | 11 | 74/84 | 15 |
| 10 | 152 | 154 | 3 | 69/71 | 2 |

At \(m=10\), suppressing degree-two paths leaves four kernel segments:
three paths from grand source node 0 to grand sink node 718, of lengths
1, 3, and 35, and one length-32 loop based at node 0. Thus the unresolved
152-terminal indicator dual has a three-cycle core. The sharp basis changes
at \(m=8,9,10\) are visible directly in the cycle ranks; no stationary
terminal kernel underlies the five tested selector margins.

That specific three-cycle routing can now be closed at the topology level.
Keep only the 152 terminal games, impose the complete exact cone and game
monotonicity at each terminal, and add precisely the worth equalities and
inequalities forced by the 154 protected source–sink path pairs. The
resulting relaxation has

\[
4864\ \text{variables},\qquad
56244\ \text{inequalities},\qquad
3708\ \text{equalities}.
\]

Its maximum lower-expectation objective is exactly zero. A rational Farkas
certificate uses 433 inequality rows and 1,067 equality rows, reconstructs
with denominator at most 100,000, and has variable objective \(-71\);
the sink constant is also \(-71\), so the original maximized objective is
zero. A common additive game at every terminal attains zero. The separate
verifier rebuilds all terminal constraints from the path report and checks
stationarity and objective over the rationals.

Therefore no change to internal games or path lengths can make this exact
\(m=10\) terminal routing positive. A negative-selector search must alter
the protected terminal incidence or its coalition labels, not merely
reoptimize the archived three-cycle kernel. This is a topology-specific
finite no-go certificate, not yet a theorem for every interlocking
indicator network.

The same exact terminal-only test also closes the \(m=6\), \(m=7\), and
\(m=9\) routings. The perspective construction below closes \(m=8\):

| \(m\) | Terminals / protected pairs | Exact dual rows (ineq./eq.) | Variable objective |
|---:|---:|---:|---:|
| 6 | 54 / 55 | 135 / 302 | \(-24\) |
| 7 | 73 / 74 | 230 / 460 | \(-34\) |
| 8 | 99 / 102 | 20,848 / 1,662 | \(-80\) |
| 9 | 152 / 163 | 375 / 912 | \(-113\) |
| 10 | 152 / 154 | 433 / 1,067 | \(-71\) |

At \(m=9\), the sole mixed-sign terminal is still affine-indicator:
\((1,-1,-1,1,-1)=2\mathbf1_9-\mathbf1_N\), so exactness linearizes its
lower expectation. Independent verification passes for all four rational
indicator-terminal certificates; a separate verifier handles the
perspective certificate.

The exceptional formulation is \(m=8\). It has two genuinely
three-level sink divergences,
\((0,-2,-1,-1,-2)\) and \((-2,0,-1,-1,-2)\). Each has 397 extreme
core-dual representations. A joint branch MILP finds a zero incumbent, but
the one-hot, binary, and hybrid big-\(M\) encodings retain upper gap
\(0.05\).

A new perspective formulation continuously convexifies all 397 branches at
each mixed sink by using branch-weighted exact-monotone copies of the
terminal game. It is an outer relaxation of every original realization.
The 29,370-variable, 323,334-inequality, 4,056-equality LP solves to zero,
and its MILP form closes at the root with reported gap zero.
Simplex-basis support reduces the certificate to a 3,230-variable rational
system of rank 3,228. Fixing the two free dual-decomposition coordinates at
\(1/2\) gives an exact rational solution. It uses 20,848 nonzero inequality
rows and 1,662 nonzero equality rows and has variable objective \(-80\),
exactly matching the sink constant. A common additive terminal realization
attains zero. An independent verifier reconstructs all 29,370 stationarity
coordinates and validates every sign and objective term, then regenerates
all 794 local representation blocks. Because every original branch realization
embeds one-hot in the perspective relaxation, this closes the \(m=8\)
terminal routing exactly.

This is a sharp weak-boundary sequence, not a counterexample.

The component gaps—not the larger global envelope gaps—are the relevant
empty-core measure. For \(m=6,7,8,9,10\) they satisfy

\[
\boxed{
g_{\rm conn}(m)=2-\frac{12}{m}.
}
\]

Thus \(m=6\) is exactly on the connected common-core boundary, while the
tested \(m>6\) families are connected empty-core obstructions. Their
normalized candidate law is

\[
\frac{t_C(m)}{g_{\rm conn}(m)}
=\frac{9m-28}{8(m-6)(16m-49)}
\sim\frac9{128m}>0.
\]

## Cross-resolution coupling

Raising several resolutions to one common grand worth creates additional
legal proper-bump edges. Unions of the separately active supports give exact
small certificates, including

\[
t_C(5{:}8)=\frac1{90}
\]

and a best player-swapped value \(2/195\). These obey a fixed-residual
three-node fork recurrence.

The separately active supports are not complete for cross-resolution work.
The full \(m=5,6\) union has 32,715 games and 155,432 edges, and compresses
to the independently verified exact archive

\[
\boxed{
\text{games/edges}=38/78,\quad
g=\frac65,\quad
t_C=\frac1{60},\quad
t_Q=\frac7{160}.
}
\]

Adding every \(m=7\) reachable state to that exact support lowers the margin
again:

\[
\boxed{
\text{games/edges}=45/160,\quad
g=\frac65,\quad
t_C=\frac1{90},\quad
t_Q=\frac7{240}.
}
\]

Both duals have the same unnormalized residual \(-4/5\); their protected
flow normalizations are 48 and 72. The next blockwise \(m=8\) calculation
was predicted to have normalization 96 and margin \(1/120\). The exact
\(m=8\) block instead has 51 games, 244 edges, and

\[
t_C=\frac2{225},\qquad t_Q=\frac7{300};
\]

its normalization is 90, so the prediction fails at an active-regime
change.

A later component audit finds active-component gap \(-4/5\) in all three
cross-resolution block archives. Their positive global gap \(6/5\) is
supplied by disconnected envelope games. These remain valid exact topology
experiments but are not connected empty-core obstructions.

## Symbolic rolling-triple theorem

The best-swap \(m=5,6,7\), \(m=6,7,8\), and compacted \(m=7,8,9\) archives
are isomorphic as coalition-labeled 43-node, 136-edge graphs. Interpolating
their games and allocations affinely in \(1/L\) gives two exact charts whose
union covers every real \(5\le L\le7\):

\[
\boxed{
v_L(N)=6+\frac4L,\quad
\operatorname{gap}(L)=2-\frac4L,\quad
t_C(L)=\frac{2}{33L}.
}
\]

`n5_parametric_rolling_triple_verify.py` checks over the entire parameter
interval:

1. game monotonicity and all 280 facets of the complete five-player exact
   cone;
2. every strict single-coordinate edge;
3. an exact parametric core primal;
4. a matching fixed rational dual; and
5. exact common-core threshold \(8\), by both an upper allocation and a
   balanced lower certificate.

The first affine chart covers \(5\le L\le6\); the second covers
\(6\le L\le7\). They meet at the exact-cone active-regime change at \(L=6\).
The second endpoint is independently verified with

\[
\text{games/edges}=43/136,\qquad
g=\frac{10}{7},\qquad
t_C=\frac{2}{231}.
\]

An exhaustive audit of 18,736 affine inequalities per chart proves that
these archived charts have maximal domains exactly \([5,6]\) and \([6,7]\).
A third chart—not extrapolation of the second—is required beyond \(L=7\).

The active selector components at the three integer supports have
common-core gaps \(-4/5,-2/3,-4/7\). Hence the symbolic theorem is correct
but its empty global common core is assembled across disconnected
components. It is not the connected weak-boundary sequence identified
above.

## Structural theorem explaining the sign

The active cross-resolution dual is a source with divergence proportional
to \((1,1,1,1,2)\) and two overlapping indicator sinks. The strengthened
theorem in `COMPLEMENTARY_FLOW_NO_GO.md` proves that any single-source
indicator fork, with arbitrary overlapping coalition blocks and weights,
has nonpositive lower-expectation objective. The common-sink reverse is
also nonpositive, and the result allows varying grand worths.

Thus these sequences can lengthen protected paths and drive the normalized
margin to zero, but their residual cannot change sign. A counterexample must
use a genuinely interlocking multi-source/multi-sink network or mixed
terminal divergences.

## Verification status

Every exact archive cited above has passed a separate verifier checking:

- rational worths and allocations;
- game monotonicity;
- all 280 exact-cone facets;
- exactly one positive worth change on each edge;
- the exact common-core gap;
- the core primal and matching rational dual;
- the coordinate-box primal and dual; and
- the non-atomic facet tax.

“Independent” means separate local verification code, not external
replication or peer review.

## Conclusion

\[
\boxed{\text{The universal coalitionally monotone selector remains open.}}
\]

The main advance is now an exact, structurally explained \(1/m\) boundary
law and a symbolic continuum theorem. Neither supplies the required
\(t_C<0\) certificate or a universal nonnegative-flow decomposition.
