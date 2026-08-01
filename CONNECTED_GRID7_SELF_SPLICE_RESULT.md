# Connected \(m=7\) self-splice search

**Date:** 2026-07-30  
**Central selector question:** open

## Starting component

The full reachable \(m=7\) support contains a standalone connected selector
component with

\[
\text{games/edges}=370/510,\qquad
g_{\rm conn}=\frac27,\qquad
t_C=\frac5{252},\qquad
t_Q=\frac3{77}.
\]

The component itself has an empty common core; its gap is not supplied by
disconnected envelope games. Its dual has 73 signed-indicator terminals
interlocked across many sources and sinks.

## First-generation connected splices

Two copies of the component were glued by identifying a sink in the first
copy with a grand-indicator source in the second. The second copy's bump
scale was \(1/10\). The exact-cone game vertex was then optimized against
the combined dual and alternated once against the true selector dual.

The best tested identification uses the sink with coalition \(25\):

\[
\boxed{
\begin{aligned}
\text{games/edges}&=739/1020,\\
g_{\rm conn}&=\frac9{161},\\
t_C&=\frac{197}{54740},\\
t_Q&=\frac{11}{2240},\\
\tau&=\frac{1149}{875840},\\
\frac{t_C}{g_{\rm conn}}&=\frac{197}{3060}.
\end{aligned}
}
\]

The archive is one connected component and passed the separate exact
verifier. Its normalized ratio \(197/3060\approx0.06438\) improves on the
base \(5/72\approx0.06944\) by about \(7.3\%\).

The next-best tested identification uses coalition \(15\):

\[
\text{games/edges}=739/1020,\quad
g_{\rm conn}=\frac{48}{805},\quad
t_C=\frac{1243}{310730},\quad
\frac{t_C}{g_{\rm conn}}=\frac{1243}{18528}.
\]

It is also one connected independently verified exact archive. Coalition
\(21\) and \(26\) identifications were worse in floating normalized margin.

The best exact dual has six signed-indicator terminals: two grand sources
and a coalition-\(10\) source feed a grand sink and sinks at coalitions
\(15\) and \(26\). This is a genuinely multi-source/multi-sink topology,
not a single fork or join covered by `COMPLEMENTARY_FLOW_NO_GO.md`.

### Scale-control audit and exact breakpoint

The original \(1/10\) scale was a search choice, not an optimum. Holding the
same sink/source identification fixed and varying the second-copy bump scale
\(s\) reveals adjacent affine regimes:

\[
\begin{array}{ll}
g(s)=\dfrac{8s}{7},&
t_C(s)=\dfrac{23000s+2625}{1368500},\\[6pt]
g(s)=\dfrac3{23}-\dfrac{4s}{7},&
t_C(s)=\dfrac{2s}{119}+\dfrac3{1564}.
\end{array}
\]

The normalized ratio decreases on the first regime and increases on the
second. Their intersection is

\[
\boxed{
s_*=\frac7{92},\qquad
g_{\rm conn}=\frac2{23},\qquad
t_C=\frac5{1564},\qquad
t_Q=\frac3{736},\qquad
\tau=\frac{11}{12512},\qquad
\frac{t_C}{g_{\rm conn}}=\frac5{136}.
}
\]

Floating solves on both sides and at the breakpoint reproduce these
formulas. More strongly, n5_splice_scale_chart_verify.py checks every
game, edge, core constraint, box constraint, and fixed primal-dual
certificate as an affine rational function of \(s\). The left certificate
has maximal archived domain

\[
\frac{37}{500}\le s\le\frac7{92},
\]

and the right certificate has maximal archived domain

\[
\frac7{92}\le s\le\frac{39}{500}.
\]

The audit covers 325,390 exact affine constraints on each chart. Thus
\(s=7/92\) is a proved local normalized optimum over the two maximal
adjacent certificate regimes, not just a floating breakpoint. The rational
breakpoint archive is exact and its one-component audit recovers the full
gap \(2/23\) in the component carrying the active dual. The control also
rules out a trivial small-bump collapse: scales
\(1/100\) and \(1/40\) plateau at ratio \(1/17\), while scale \(1/5\) is
infeasible.

Recursing from the tuned archive at the same outer scale gives the exact
1,477-game, 2,040-edge family:

\[
g_{\rm conn}=\frac2{23},\qquad
t_C=\frac{467}{207368},\qquad
t_Q=\frac{153}{35443},\qquad
\tau=\frac{28687}{13893656},\qquad
\frac{t_C}{g_{\rm conn}}=\frac{467}{18032}.
\]

Its normalized value \(467/18032\approx0.02590\) is \(29.6\%\) below
the tuned first splice and only about \(11.4\%\) above the full-\(m=9\)
record \(53/2280\). The sparse rational generator and independent verifier
both pass, as does the one-component gap audit. Its dual has 49 active flow
rows and eight signed-indicator terminals.

A third tuned module gives the new exact normalized record:

\[
\boxed{
\begin{aligned}
\text{games/edges}&=2215/3060,\\
g_{\rm conn}&=\frac2{23},\\
t_C&=\frac{237}{135424},\\
t_Q&=\frac{2567}{368184},\\
\tau&=\frac{2675}{512256},\\
\frac{t_C}{g_{\rm conn}}&=\frac{237}{11776}.
\end{aligned}
}
\]

The ratio \(237/11776\approx0.02013\) is \(22.3\%\) below tuned level two
and \(13.4\%\) below the former exact normalized record \(53/2280\) from
the full \(m=9\) component. The sparse exact archive, separate verifier,
and one-component audit all pass. Its dual has scale 64, 64 active flow
rows, ten signed-indicator terminals, and a new coalition-\(15\) sink at
node 1585. The simple residual recurrence inferred from tuned levels one
and two changes basis here, but the normalized margin falls further.

A fourth tuned module is also exact:

\[
\boxed{
\begin{aligned}
\text{games/edges}&=2953/4080,\\
g_{\rm conn}&=\frac2{23},\\
t_C&=\frac{481}{334328},\\
t_Q&=\frac{13219}{4418208},\\
\tau&=\frac{542137}{349038432},\\
\frac{t_C}{g_{\rm conn}}&=\frac{481}{29072}.
\end{aligned}
}
\]

Its normalized ratio \(481/29072\approx0.01655\) is \(17.8\%\) below tuned
level three. The sparse certificate initially required a joint optimal-face
repair because independent coordinate rounding did not preserve its
degenerate active equalities. The repaired exact archive and the independent
verifier both pass; its selector graph is one component with the full
\(2/23\) gap. Its dual has scale 79, 79 active flow rows, 12 terminals, and
a new coalition-\(15\) sink at node 2323.

A fifth tuned module changes basis again:

\[
\boxed{
\begin{aligned}
\text{games/edges}&=3691/5100,\\
g_{\rm conn}&=\frac2{23},\\
t_C&=\frac1{828},\\
t_Q&=\frac{1435}{622633},\\
\tau&=\frac{24589}{22414788},\\
\frac{t_C}{g_{\rm conn}}&=\frac1{72}.
\end{aligned}
}
\]

The ratio \(1/72\approx0.01389\) is \(16.1\%\) below tuned level four and
\(40.3\%\) below the former full-\(m=9\) record \(53/2280\). The sparse
exact archive, independent verifier, and one-component audit all pass. Its
dual has scale 90, 90 active flow rows, two active core rows, and 13
terminals. The final proper-coalition terminal pivots from a
coalition-\(15\) sink to a coalition-\(5\) sink at node 3061.

Although this 13-terminal dual fails both orientations of the general
block-transport audit, it has an exact alternating-ladder decomposition.
Six grand sources split the complementary player blocks \(5\) and \(26\)
through five grand sinks and the two proper endpoint sinks. Protected-path
invariance makes the relevant block worth nonincreasing along the ladder;
the endpoint objective is bounded by the first source's complementary-pair
balancedness inequality. The exact ladder audit gives unnormalized
objective \(-5/46\) and safe slack \(5/46\). Thus the selected fifth-level
dual is structurally safe even though its value falsifies the earlier
recurrence.

A sixth tuned module stays in the alternating-ladder regime and gives the
current exact normalized record:

\[
\boxed{
\begin{aligned}
\text{games/edges}&=4429/6120,\\
g_{\rm conn}&=\frac2{23},\\
t_C&=\frac{473}{444360},\\
t_Q&=\frac{1691}{414736},\\
\tau&=\frac{18743}{6221040},\\
\frac{t_C}{g_{\rm conn}}&=\frac{473}{38640}.
\end{aligned}
}
\]

The ratio \(473/38640\approx0.0122412\) is \(11.86\%\) below level five
and \(47.34\%\) below the former full-\(m=9\) record. The sparse exact
archive, independent verifier, and one-component audit all pass. Its dual
has scale 105, 105 active flow rows, two active core rows, and 15
signed-indicator terminals: seven grand sources alternate through six grand
sinks to endpoint sinks at coalitions 26 and 5.

The exact alternating-ladder audit gives unnormalized objective
\(-473/4232\) and matching safe slack \(473/4232\). Independently, the
terminal-only relaxation has exact optimum zero; its rational certificate
uses eight inequalities and 27 equalities. Thus this sixth selected basis
is closed both by the ladder theorem and by a finite exact
terminal-topology certificate. Lengthening this unchanged ladder cannot
cross zero.

The first four exact tuned levels satisfy

\[
\boxed{
\frac{t_k}{g_{\rm conn}}
=\frac{7k+453}{368(15k+19)},\qquad
t_k=\frac{7k+453}{4232(15k+19)}
}
\qquad (k=1,2,3,4).
\]

The dual scale law \(15k+19\) and terminal-count law \(2k+4\) also hold at
those four levels. Level five falsifies every displayed continuation:
the old formula predicts \(61/4324\), whereas the exact answer is \(1/72\);
the scale is 90 rather than 94, and there are 13 rather than 14 terminals.
The conjectured positive limit \(7/5520\) therefore does not retire this
route. Level six confirms that the new two-core-row regime can keep lowering
the margin, but its selected topology is an exactly safe alternating ladder.
Future progress on this branch therefore requires another basis or terminal
topology change, not merely another copy of the same ladder.

## Second-generation exact splice

Gluing two copies of the best 739-game archive at its coalition-\(15\) sink
produces one connected exact family:

\[
\boxed{
\begin{aligned}
\text{games/edges}&=1477/2040,\\
g_{\rm conn}&=\frac9{161},\\
t_C&=\frac{401}{157780},\\
t_Q&=\frac{272}{53935},\\
\tau&=\frac{5289}{2114252},\\
\frac{t_C}{g_{\rm conn}}&=\frac{401}{8820}.
\end{aligned}
}
\]

The component audit finds exactly one weak component containing all 1,477
games, all 2,040 edges, the active dual flow, and the full gap \(9/161\).
The separate archive verifier independently recomputes the family optimum
and passes.
The normalized ratio \(401/8820\approx0.04546\) is about \(29.4\%\) below
the first-generation ratio and \(34.5\%\) below the base component.

The exact dual has scale \(49\), 49 active flow rows, three active core
rows, and eight signed-indicator terminals. A path from node 92 to node 836
crosses the splice through the identified edge \(109\to823\), so the
certificate is not merely the disjoint sum of the two parent certificates.
It has three grand sources, a coalition-\(10\) source, two grand sinks, a
coalition-\(26\) sink, and a coalition-\(15\) sink.

## Third-generation exact splice and basis change

Attaching a third 739-game module at the new coalition-\(15\) sink gives a
one-component 2,215-game, 3,060-edge exact optimum:

\[
g_{\rm conn}=\frac{97}{1610},\qquad
t_C=\frac{51}{25760},\qquad
t_Q=\frac{1201}{280140},\qquad
\tau=\frac{5171}{2241120},\qquad
\frac{t_C}{g_{\rm conn}}=\frac{51}{1552}.
\]

The normalized value \(51/1552\approx0.03286\) is about \(27.7\%\) below
generation two and \(52.7\%\) below the base component.

This out-of-sample point falsifies the first affine-recurrence guess from
generations one and two. That guess predicted \(g=9/161\) and normalized
margin \(121/3456\); instead, the exact-cone optimum moves to a new basis
with the larger gap \(97/1610\) and a smaller ratio. The sparse exact
archive passes all 280 exact-cone facets, all 3,060 exact one-coordinate
edge checks, exact core and box primal-dual stationarity, and the separate
sparse archive verifier. The selector graph is one component, so its exact
\(97/1610\) gap and the active dual are component-aligned.

The common-core certificate itself changes regime. Generations one and two
use the complementary envelope pair \(3,28\). Generation three uses

\[
v(12)=\frac{33}{115},\qquad
v(19)=\frac{249}{322},\qquad
v(12)+v(19)-v(N)=\frac{97}{1610}.
\]

Thus the improvement is accompanied by a genuine complementary-envelope
pivot, not merely a longer denominator in the old basis.

## Interpretation

This is the first recursive coupling in the continuation that simultaneously

1. keeps the empty-common-core witnesses in the active connected component;
2. lowers the normalized selector margin;
3. retains positive non-atomic facet tax; and
4. avoids collapse to a single indicator star.

Every verified margin remains strictly positive. The construction therefore
advances the connected frontier but does not settle the universal-selector
question.

## Reproducible exact outputs

- n5_sparse_exact_certificate.py builds the sparse rational primal-dual
  archives, including the 4,429-game sixth tuned splice.
- n5_verify_sparse_exact_archive.py independently verifies their game,
  edge, common-core, core-margin, and box-margin certificates.
- n5_splice_scale_chart_verify.py proves the two adjacent affine scale
  charts and their maximal archived domains.
- n5_indicator_transport_audit.py proves that the 6-, 8-, and 10-terminal
  splice duals do not reduce to the block-transport theorem.
- n5_alternating_ladder_audit.py verifies the exact fifth- and sixth-level
  alternating-ladder decompositions.
- n5_terminal_topology_relaxation.py and
  n5_verify_terminal_topology_certificate.py give a second exact no-go
  certificate for the sixth-level terminal routing.
