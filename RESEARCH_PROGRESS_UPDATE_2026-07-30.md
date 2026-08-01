# Research progress update: universal coalitionally monotone core selection

**Updated:** 2026-07-30  
**Project:** `economics-research/game-theory/exact-game-monotone-selection`  
**Central status:** **open**

No finite exact-game family with compatibility margin \(t_C\le0\) has been
verified, and no universal existence theorem has been proved. The work has,
however, produced a new exact empty-common-core frontier, eliminated several
plausible topology classes, and isolated the remaining obstruction much more
sharply.

## Basic explanation

A cooperative game gives every group of players \(S\) a value \(v(S)\). A
core allocation divides the total value while giving every group at least
its own value. An exact game is one whose core can make every group’s
constraint tight somewhere.

The desired universal rule would choose one core allocation from every exact
game. If only one coalition becomes more valuable, none of that coalition’s
members may lose.

For any finite network of games and one-coalition increases, this is an exact
linear-programming question. Its max-min compatibility margin is \(t_C\):

| Margin | Meaning |
|---:|---|
| \(t_C>0\) | This finite family has a strictly compatible choice. |
| \(t_C=0\) | Weak compatibility survives exactly on the boundary. |
| \(t_C<0\) | No compatible choice exists; this disproves a universal selector. |

Tiny positive values are not counterexamples. The required central result is
a rationally verified \(t_C\le0\), or a theorem proving \(t_C\ge0\) for
every finite family.

There is one further interpretation rule. Monotonicity constraints decompose
over weakly connected components of the game graph. A positive common-core
gap computed across disconnected components can coexist with a common core
inside every selector component. The component carrying the active dual must
therefore be audited separately; global gap-to-margin ratios are not used
when the envelope games are disconnected from it.

## Current exact frontier

The strongest independently verified connected empty-common-core family now
has 219 games and 218 legal bump edges:

\[
g=\frac{37}{140{,}000{,}000}>0,
\qquad
t_C=\frac{1403}{2{,}240{,}000{,}000}
\approx6.2634\times10^{-7}>0.
\]

It was built by attaching both orientations of a verified 91-edge
complement path to the exact 37-game non-large-core family. The separate
archive verifier reports:

- all 219 games exact;
- all 218 edges legal one-coordinate increases;
- empty common core with the displayed exact gap; and
- an exact primal-dual monotonicity optimum with the displayed margin.

The one-sided predecessor is also independently verified:

\[
\begin{aligned}
\text{games/edges}&=128/127,\\
g&=\frac9{350{,}000{,}000},\\
t_C&=\frac{43{,}849}{50{,}400{,}000{,}000}
\approx8.7002\times10^{-7}.
\end{aligned}
\]

These margins are much smaller than the earlier transverse certificates:

| Exact family | Games / edges | Empty-core gap | Exact margin |
|---|---:|---:|---:|
| Correct transverse stage 2 | 310 / 917 | \(5329/44{,}800{,}000\) | \(69/946{,}400\) |
| Correct transverse stage 3 | 356 / 1,109 | \(5409/44{,}800{,}000\) | \(69/1{,}388{,}800\) |
| Empty pair/triple–nonlarge hybrid | 316 / 923 | \(1/10{,}000\) | \(100{,}993/848{,}000{,}000\) |
| Empty/taxed recursive hybrid | 322 / 929 | \(1/10{,}000\) | \(101{,}003/856{,}000{,}000\) |
| One-sided complement path | 128 / 127 | \(9/350{,}000{,}000\) | \(43{,}849/50{,}400{,}000{,}000\) |
| Two-sided complement paths | 219 / 218 | \(37/140{,}000{,}000\) | \(1403/2{,}240{,}000{,}000\) |

The 219-game result is the main new exact research output. Its coordinate-box
margin is \(1009/432{,}000\), so nearly all of the box slack is consumed by
genuinely non-atomic exact-core facets.

## Connected structural frontier: cross-resolution overlap

The raw \(219\)-game margin above is numerically smaller, but a new
cross-resolution construction supplies a second meaningful normalization:
the component carrying its active selector dual has common-core gap \(1\),
and its worths are not shrunk toward a scaling degeneracy.

Canonical intrinsic paths at resolutions \(m=5,6,7\) were joined to the
previous \(m=2,3,4\) family on the common grand-worth section \(v(N)=7\).
Every retained state was checked against all 280 exact-cone facets. The
resulting exact sequence is:

| Resolutions | Games / edges | Exact margin | Box margin | Facet tax |
|---|---:|---:|---:|---:|
| \(2,\ldots,5\) | 402 / 627 | \(1/36\) | \(3/98\) | \(5/1764\) |
| \(2,\ldots,6\) | 520 / 763 | \(1/39\) | \(1/34\) | \(5/1326\) |
| \(2,\ldots,7\) | 658 / 917 | \(1/42\) | \(1/36\) | \(1/252\) |

The active selector components in the three full families have respectively
297/496, 415/632, and 553/786 games/edges, and each has exact common-core
gap \(1\). Thus the full-family sequence survives the component audit.

The \(m\le7\) archive can be compressed to 17 games and 57 induced edges
without changing the *global* gap or margin. Adding its best player-swapped
image produces a separately verified 21-game, 90-edge family with

\[
g=1,\qquad t_C=\frac1{54},\qquad
t_Q=\frac1{24},\qquad \tau=\frac5{216}.
\]

However, compression disconnects the envelope witnesses: the 19-game active
component of the swapped archive has gap \(-1\). The improved \(1/54\)
margin is therefore a topology side result, not a stronger connected
normalized frontier. The uncompressed 553-game active component with
\(g_{\rm conn}=1\) and \(t_C=1/42\) is the accepted connected certificate.

Its exact dual has 54 active monotonicity rows, all weighted \(1/54\).
After removing this normalization, the source lower expectation is \(7\)
and the two sink expectations are \(-4\) each. Thus the residual is the
fixed integer \(-1\); resolution and the productive swap lengthen the
protected paths but do not change its sign.

All 60 player-permutation images of the \(m\le6\) compact family, all 79
admissible modular alignments, and affine-modular self-compositions were
screened. They plateau at \(1/48\), the best single-swap value for that
resolution. This closes simple symmetry saturation on the common
grand-\(7\) section.

## Stronger floating leads, with the necessary caveat

Attaching proper 218-game bridge gadgets to the two-sided family gives
positive floating margins below \(10^{-7}\). The best current asymmetric
chain is:

| Stage | Floating margin |
|---:|---:|
| 1 | \(1.2027\times10^{-7}\) |
| 2 | \(6.9757\times10^{-8}\) |
| 3 | \(4.9125\times10^{-8}\) |
| 4 | \(4.3060\times10^{-8}\) |

These are discovery results, not exact certificates. Their active dual has a
constant positive numerator and an increasing normalization denominator.
Moreover, arbitrarily small bump scales can make a strict margin small
without approaching weak incompatibility in a meaningful normalized sense.
They therefore must not be described as “almost counterexamples” solely
because of their decimal size.

The first 436-game exactification attempt reconstructed exact games and the
box optimum but failed to reconstruct the exact core optimum. A second run
with a billion-scale rational denominator failed at the same core
reconstruction step. The sub-\(10^{-7}\) chain therefore remains
provisional and is not part of the certified frontier.

## Structural conclusions from the latest search

### 1. Arbitrary indicator forks and joins are structurally nonpositive

The recurring active dual uses complementary coalitions, such as 12 and 19.
Appending verified bridges changes the dual normalization:

\[
t_C=\frac{\text{fixed positive residual}}{\text{growing positive weight}}.
\]

Two-ended path closures, all 12 internal player orientations, a two-bridge
algebraically closed ring, and four additional aligned stages all remained
positive.

There is also a proof-level reason. A protected flow for player \(i\) can use
only edges whose changed coalition contains \(i\). Along such an edge, every
coalition excluding \(i\) has unchanged worth.

This gives a general no-go theorem. If one source routes any weighted
collection of coalition indicators to separate sinks, exactness moves each
sink's complementary worth back to the source. Evaluating the source lower
expectation at any one of its core points then proves the total dual
objective is nonpositive. The reverse common-sink join follows by the dual
argument. No balanced-cover hypothesis is needed. Complementary pairs,
overlapping blocks, player-partition forks, co-singleton stars, the active
\(2,5,24\) mechanism in the 322-game hybrid, and the new cross-resolution
fork are all special cases.

The strengthened proof is recorded in `COMPLEMENTARY_FLOW_NO_GO.md`. It
explains the constant-numerator recurrences and retires every indicator-flow
star, a substantially larger topology class than the complementary chains
alone.

### 2. Direct active/envelope alignment is structurally blocked

For the 322-game hybrid, the margin dual is supported at coalitions
\(2,5,24\), a genuine three-part partition, while the empty common core is
witnessed by the five co-singletons. For the asymmetric tiny-margin chain,
the active pair is at different nodes from the two complementary envelope
maxima.

LP tests that force the active nodes themselves to attain the empty-envelope
values are infeasible. Exact coordinate-closure and two-ended bridge tests
either become infeasible or move to a different positive dual.

### 3. Balanced stars alone are positive

The exact five-leaf co-singleton star has

\[
g=\frac1{80},
\qquad
t_C=\frac1{400}.
\]

A complete one-step star containing every proper coalition optimized to the
floating formula \(t_C=\delta/4\); adding the grand-coalition child gave
\(t_C=\delta/5\). These are strong computational indications that a
single branching layer is not enough. A counterexample needs branch-and-
rejoin geometry.

### 4. Balanced non-atomic gluing is possible but still positive

The independently reverified 55-game legal grid has a taxed four-player
active port. Grafting one copy onto a co-singleton-star leaf lowered the
margin from \(0.0025\) to \(0.0020858333\) while preserving an empty-core gap
near \(0.012495\). Five leaf copies and a stabilizer-orbit copy did not
accumulate the tax: the dual selected one copy and ignored the others.

The five-copy family has:

\[
t_C\approx0.0020858333,\quad
t_Q\approx0.0020875,\quad
\tau=t_Q-t_C\approx1.67\times10^{-6}.
\]

This is the correct non-complementary mechanism, but it is far too weak in
its current uncoupled form.

### 5. The crown and aggregate branches are closed for the tested templates

The crown aggregate bridge remained positive throughout its feasible scale
range; near the boundary its margins rose from about \(4.03\) to
\(4.15\times10^{-7}\), and larger scales were infeasible. Complete and
partial cube closures also remained positive.

### 6. Reproduction of the compact exact overlapping-flow backbone

The 304-game legal cross-resolution family was independently reduced to the
active dual nodes plus two singleton-envelope nodes. Taking every induced
legal edge reproduces, up to node ordering, the existing project archive
`n5_intrinsic_m234_grand7_compact_exact.json`:

\[
\begin{aligned}
\text{games/edges}&=77/124,\\
g&=1,\\
t_C&=\frac4{133},\\
t_Q&=\frac3{94},\\
\tau&=\frac{23}{12502}.
\end{aligned}
\]

All quantities are identical to the 304-game parent. The independent
reconstruction and verifier both pass. The compact archive’s optimal dual
uses indicator divergences, correcting the initial impression that mixed
coefficients were essential. What remains structurally interesting is the
overlap of several indicator joins/forks at shared nodes; the network is not
one balanced fork or join covered by the new no-go theorem. This 77-game
family remains a useful topology laboratory, even though its raw margin is
not the numerical frontier.

The later selector-component audit finds three weak components. The
75-game component carrying the active dual has common-core gap \(-1/3\);
the global gap \(1\) comes from the two isolated envelope games. Recursive
self-splices can make the displayed margin extremely small, but their active
components also retain nonempty common cores, and their sparse duals collapse
to balanced or complementary forks. They are therefore retired as
counterexample routes. Details are in `CONNECTED_COMPONENT_AUDIT.md`.

### 7. A sharp lower-grand path sequence and its exact-cone barrier

The grand-\(7\) paths were not minimax in grand worth. Exhaustive reachable
state searches at resolutions \(m=5,6,7\) give the sharp discrete caps

\[
G_m=6+\frac4m:
\qquad
G_5=\frac{34}{5},\quad
G_6=\frac{20}{3},\quad
G_7=\frac{46}{7}.
\]

At each resolution, the immediately lower lattice cap fails at depth 12.
Every exit from the reachable component has exact root at least \(G_m\).
The same eight blocking states and the same exact-facet types
9, 10, 11, 13, and 14 recur in all three cases. This is an exact finite
graph calculation with rational facet evaluations, not a floating
infeasibility judgment.

The independently verified canonical-path archives are:

| Resolution | Games / edges | Gap | Exact margin | Box margin |
|---:|---:|---:|---:|---:|
| 5 | 101 / 100 | \(6/5\) | \(7/90\) | \(2/15\) |
| 6 | 121 / 120 | \(4/3\) | \(11/162\) | \(11/96\) |
| 7 | 141 / 140 | \(10/7\) | \(23/378\) | \(9/91\) |

The common-core gap follows \(2-4/m\). More surprisingly, the exact core
margins fit

\[
t_C(m)=\frac{m+16}{54m}
\]

in all three cases. The \(m=5\) and \(m=6\) optimal duals have identical
weights and support shape, and \(m=7\) agrees out of sample. This path
sequence therefore appears to approach the positive value \(1/54\), not
zero, even though its grand cap approaches the original obstruction's
grand worth \(6\).

The best player-swap pairs at \(m=5,6,7\) are also exact:

\[
\frac{42}{575},\qquad \frac{22}{345},\qquad \frac2{35}.
\]

They fit

\[
t_{\rm pair}(m)=\frac{2(m+16)}{115m},
\]

with positive limit \(2/115\). At \(m=5\), three greedily coupled
permutation images lower the margin further to \(19/290\) on 91 games and
102 edges; two additional images do not improve it.

### 8. Complete reachable components reveal an asymptotic boundary law

The canonical paths omit almost all exact states available below the sharp
grand cap. Sparse LPs now include every reachable state and every induced
legal edge. Their active dual supports were then compressed, rationally
reconstructed, and checked by the independent verifier:

| \(m\) | Full games / edges | Active component | Component gap | \(t_C\) | \(t_Q\) | Facet tax |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 13,755 / 64,916 | 195 / 291 | \(-2/15\) | \(17/620\) | \(7/135\) | \(409/16740\) |
| 6 | 18,963 / 90,448 | 394 / 584 | \(0\) | \(13/564\) | \(2/45\) | \(181/8460\) |
| 7 | 24,171 / 115,980 | 370 / 510 | \(2/7\) | \(5/252\) | \(3/77\) | \(53/2772\) |
| 8 | 29,379 / 141,512 | 494 / 679 | \(1/2\) | \(11/632\) | \(3/76\) | \(265/12008\) |
| 9 | 34,587 / 167,044 | 897 / 1,409 | \(2/3\) | \(53/3420\) | \(13/396\) | \(163/9405\) |

The full-family optima and their exact active components agree. All five
exact points satisfy

\[
t_C(m)=\frac{9m-28}{4m(16m-49)}
\]

These are observed core-margin formulas, not yet all-\(m\) theorems. They
are supported by five independently verified exact resolutions, by the
exact linear full-family counts

\[
\#V_m=5208m-12285,\qquad
\#E_m=25532m-62744
\]

for \(m=5,\ldots,9\), and by active core-dual normalization
\(8(16m-49)\). The \(m=9\) dual again has \(12m-40=68\) active core
rows. Its selected sparse basis has 123 nonzero-divergence nodes rather than
the \(22m-79=119\) continuation of the earlier bases, so the margin law
survives a second active-topology change. The formula implies

\[
t_C(m)\sim\frac{9}{64m}>0.
\]

Thus completing the entire reachable component destroys the positive floor
of the canonical path and approaches weak compatibility at rate \(1/m\),
but it still does not cross zero at any tested finite resolution.

The component audit reveals the exact connected phase transition that the
earlier global envelope gap obscured. For \(m=6,7,8,9\),

\[
g_{\rm conn}(m)=2-\frac{12}{m}.
\]

The active component is on the common-core boundary at \(m=6\) and has an
empty common core for every tested \(m>6\). Consequently the meaningful
connected normalized law is

\[
\frac{t_C(m)}{g_{\rm conn}(m)}
=\frac{9m-28}{8(m-6)(16m-49)}
\sim\frac9{128m}>0.
\]

Standalone \(m=7,8,9\) active-component archives have been rationally
rebuilt and independently verified. The \(m=6\) boundary archive is also
exactly rebuilt; the separate obstruction verifier intentionally requires a
strictly positive gap, so its zero-gap status is checked by the exact
component audit rather than reported as an “empty-core PASS.”

The first three box margins suggested
\((m+2)/(3m(m+4))\), but the independently verified \(m=8\) value is
\(3/76\), not \(5/144\). Thus the coordinate-box LP has an active-set
breakpoint at \(m=8\). The core-margin recurrence survives this out-of-sample
test; the box recurrence does not.

The \(m=5\) support was screened against all 119 nonidentity player
permutations, and three greedy symmetry-closure rounds. None lowered
\(17/620\). This indicates that its gain is caused by the full intrinsic
state geometry rather than an omitted simple symmetry copy.

The large LP originally exhausted memory when represented as dense Python
rows. It was rebuilt as a sparse COO/CSR model and first checked by
reproducing the known \(7/90\) canonical-path result. The exact supports are
therefore downstream of a separately tested implementation change, not of a
different mathematical relaxation.

### 9. Cross-resolution coupling is exact but envelope-disconnected

The sharp-cap families have different grand worths. Raising their grand
coordinates to the common value \(34/5\) preserves exactness and makes
cross-resolution proper-bump edges legal. Taking every induced edge among
the individual exact supports gives:

| Resolutions at \(v(N)=34/5\) | Exact support | Gap | \(t_C\) | \(t_Q\) |
|---|---:|---:|---:|---:|
| \(5,6\) | 35 / 60 | \(6/5\) | \(2/105\) | \(1/20\) |
| \(5,6,7\) | 40 / 106 | \(6/5\) | \(4/285\) | \(9/260\) |
| \(5,6,7\) plus best swap | 43 / 136 | \(6/5\) | \(2/165\) | \(7/220\) |
| \(5,6,7,8\) | 45 / 161 | \(6/5\) | \(1/90\) | \(1/35\) |
| \(5,6,7,8\) plus best swap | 58 / 252 | \(6/5\) | \(2/195\) | \(7/370\) |

All five rows are independently verified exact primal-dual certificates.
The later component audit changes their interpretation: the active selector
components have common-core gaps \(-4/5\) for the \(m=5,6,7\) best swap and
\(-2/3\) for the \(m=6,7,8\) best swap. Their positive displayed gaps come
from envelope games in other components. Thus these are exact topology
experiments, but each active component admits a constant common-core
selection and is not a connected obstruction.

For the unswapped unions through maximum resolution \(M=6,7,8\), the dual has
three core rows, a fixed unnormalized residual \(4/5\), and respectively
42, 57, and 72 equally weighted monotonicity rows. This gives the candidate
recurrence

\[
t_{5:M}=\frac{4}{5(15M-48)}
\qquad(M\ge6),
\]

which predicted \(t_{5:8}=1/90\) before the exact \(M=8\) solve and passed
that out-of-sample test. The best swaps at \(M=7,8\) preserve the same
residual and raise the normalization to 66 and 78. All 60 player images at
\(M=7\) screened together remain at \(2/165\) in the floating sparse LP, so
simple symmetry saturation again plateaus.

The rolling \(m=7,8,9\) union at its own common grand worth \(46/7\) has an
independently verified best-swap support of 43 games and 136 edges with

\[
g=\frac{10}{7},\qquad
t_C=\frac{2}{231},\qquad
t_Q=\frac1{35},\qquad
\tau=\frac{23}{1155}.
\]

The first dual extraction contained one redundant game and two redundant
edges; a second sparse-support solve removed them without changing either
optimum.

At the global-archive level this cross-resolution construction has a smaller
margin than joining the resolutions at their own grand worths. That
varying-grand union has an exact margin
\(2/105\), but its dual is only the trivial aggregate bottleneck: one grand
increase of \(2/21\) divided among five players. It has zero non-atomic
facet tax and cannot cross zero.

The table above couples the nodes active in the separate single-resolution
LPs. A complete-state check at \(m=5,6\) includes all 32,715 distinct games
and 155,432 induced edges. It lowers the margin further and compresses to
the independently verified exact certificate

\[
\text{games/edges}=38/78,\qquad
g=\frac65,\qquad
t_C=\frac1{60},\qquad
t_Q=\frac7{160},\qquad
\tau=\frac{13}{480}.
\]

Its dual is still a three-node overlapping fork with residual \(-4/5\), but
its protected paths have normalization 48 rather than 42. Consequently the
\(4/[5(15M-48)]\) formula is a theorem candidate for the active-support
unions only. Its 18-game active component has common-core gap \(-4/5\), so
the global gap \(6/5\) does not normalize the selector obstruction.

Adding every \(m=7\) reachable state to this exact 38-game seed lowers the
margin again and gives the independently verified blockwise certificate

\[
\text{games/edges}=45/160,\qquad
g=\frac65,\qquad
t_C=\frac1{90},\qquad
t_Q=\frac7{240}.
\]

The fork residual remains \(-4/5\), while the normalization rises from 48
to 72. Adding every \(m=8\) state produces a separately verified 51-game,
244-edge exact support with

\[
g=\frac65,\qquad
t_C=\frac2{225},\qquad
t_Q=\frac7{300},\qquad
\tau=\frac{13}{900}.
\]

This falsifies the predicted \(1/120\): the normalization is 90, not 96.
All three blockwise active components have gap \(-4/5\). These results
diagnose an active-regime change but do not form a connected empty-core
frontier.

### 10. A two-chart symbolic rolling-triple boundary theorem

The best-swap \(m=5,6,7\) and \(m=6,7,8\) archives have the same
coalition-labeled 43-node, 136-edge graph. After the verified node
isomorphism, interpolate every worth and archived allocation affinely in
\(1/L\). For every real \(5\le L\le6\), this gives:

\[
v_L(N)=6+\frac4L,\qquad
\text{common-core gap}=2-\frac4L,\qquad
t_C(L)=\frac{2}{33L}.
\]

The script `n5_parametric_rolling_triple_verify.py` checks symbolically:

- monotonicity and all 280 exact-cone facets for every game over the whole
  parameter interval;
- strict positivity and one-coordinate legality of every edge;
- a parametric core allocation attaining \(2/(33L)\);
- the fixed rational dual proving the matching upper bound; and
- a common-budget upper point and balanced lower certificate proving the
  exact threshold \(8\).

The independently verified \(m=7,8,9\) support has the same labeled graph
after its redundant node is removed. A second affine chart passes the same
symbolic checks for every real \(6\le L\le7\). The two charts meet at
\(L=6\), although their node identifications differ because an exact-cone
facet changes active regime there. Together they prove the displayed
formulas, including

\[
\boxed{t_C(L)=\frac{2}{33L}},
\]

for every real \(5\le L\le7\). This is a proof over a continuum of exact
families, not a finite list of floating samples. An exhaustive audit of
18,736 affine inequalities per chart proves that their maximal archived
domains are exactly \([5,6]\) and \([6,7]\). The endpoints are genuine
core or exact-facet regime changes. Extending the theorem toward
\(L\to\infty\) therefore requires additional charts, not blind
extrapolation. The proof and reproducible outputs are summarized in
`PARAMETRIC_ROLLING_TRIPLE_RESULT.md`.

The component audit adds an important limitation. The selector component
gaps on the three integer rolling supports are respectively
\(-4/5,-2/3,-4/7\). The symbolic theorem remains correct, but its empty
global common core is assembled across disconnected components. It is a
parametric exact-LP side theorem, not a connected weak-obstruction theorem.

### 11. Connected \(m=7\) self-coupling improves the local frontier

The standalone \(m=7\) active component has 370 games, 510 edges, connected
gap \(2/7\), and margin \(5/252\). Unlike the compressed cross-resolution
families, its empty-common-core witnesses lie inside the selector component.

Gluing two copies at an active coalition-\(25\) sink and a grand-indicator
source, with second-copy bump scale \(1/10\), gives an independently verified
one-component exact archive:

\[
\begin{aligned}
\text{games/edges}&=739/1020,\\
g_{\rm conn}&=\frac9{161},\\
t_C&=\frac{197}{54740},\\
t_Q&=\frac{11}{2240},\\
\frac{t_C}{g_{\rm conn}}&=\frac{197}{3060}.
\end{aligned}
\]

This improves the base normalized ratio from \(5/72\) to
\(197/3060\), about \(7.3\%\). Its optimal dual has several sources and
sinks and is not one indicator star. Gluing another copy at the new
coalition-\(15\) sink gives the exact second-generation record

\[
\begin{aligned}
\text{games/edges}&=1477/2040,\\
g_{\rm conn}&=\frac9{161},\\
t_C&=\frac{401}{157780},\\
t_Q&=\frac{272}{53935},\\
\tau&=\frac{5289}{2114252},\\
\frac{t_C}{g_{\rm conn}}&=\frac{401}{8820}.
\end{aligned}
\]

Its selector graph is one component, so the entire gap and active dual flow
coincide. The ratio is about \(29.4\%\) below the first-generation splice
and \(34.5\%\) below the base \(m=7\) component. Its exact dual has eight
signed-indicator terminals and contains a path crossing the identified
splice node, rather than decomposing into two independent parent
certificates. The construction and verification paths are recorded in
`CONNECTED_GRID7_SELF_SPLICE_RESULT.md`.

A scale-control audit of the first splice finds that \(1/10\) is not
canonical. Two adjacent affine regimes meet at

\[
s_*=\frac7{92},\qquad
g_{\rm conn}=\frac2{23},\qquad
t_C=\frac5{1564},\qquad
t_Q=\frac3{736},\qquad
\tau=\frac{11}{12512},\qquad
\frac{t_C}{g_{\rm conn}}=\frac5{136}.
\]

The normalized ratio falls on the left chart and rises on the right, making
\(7/92\) the local scale optimum for this fixed identification. Smaller
scales \(1/100\) and \(1/40\) plateau at ratio \(1/17\), so the improvement
is not a trivial vanishing-bump effect. Scale \(1/5\) is exact-cone
infeasible. The breakpoint archive is exact, one connected component, and
has the full \(2/23\) gap in the active component.

Recursing from the tuned module at the same scale gives the exact
second-level family

\[
g_{\rm conn}=\frac2{23},\qquad
t_C=\frac{467}{207368},\qquad
t_Q=\frac{153}{35443},\qquad
\tau=\frac{28687}{13893656},\qquad
\frac{t_C}{g_{\rm conn}}=\frac{467}{18032}.
\]

The ratio \(467/18032\approx0.02590\) is \(29.6\%\) below the tuned first
splice and close to the current full-\(m=9\) normalized record. Its sparse
exact primal-dual archive, separate verifier, and one-component audit all
pass.

A third tuned module gives the new exact normalized frontier:

\[
g_{\rm conn}=\frac2{23},\qquad
t_C=\frac{237}{135424},\qquad
t_Q=\frac{2567}{368184},\qquad
\tau=\frac{2675}{512256},\qquad
\frac{t_C}{g_{\rm conn}}=\frac{237}{11776}.
\]

The exact ratio \(237/11776\approx0.02013\) is \(22.3\%\) below tuned level
two and \(13.4\%\) below the former full-\(m=9\) record \(53/2280\). Its
sparse archive, separate verifier, and one-component audit all pass. The
dual has scale 64, 64 active flow rows, and ten terminals. The simple
residual recurrence inferred from the first two tuned levels changes basis,
but the connected mechanism improves further.

A fourth tuned module is exact:

\[
g_{\rm conn}=\frac2{23},\qquad
t_C=\frac{481}{334328},\qquad
t_Q=\frac{13219}{4418208},\qquad
\tau=\frac{542137}{349038432},\qquad
\frac{t_C}{g_{\rm conn}}=\frac{481}{29072}.
\]

The ratio \(481/29072\approx0.01655\) is \(17.8\%\) below tuned level
three. Its repaired sparse primal-dual archive, independent verifier, and
one-component audit pass. Its dual has scale 79, 79 active flow rows, and
12 terminals.

A fifth tuned module gives the exact family

\[
\text{games/edges}=3691/5100,\qquad
g_{\rm conn}=\frac2{23},\qquad
t_C=\frac1{828},\qquad
t_Q=\frac{1435}{622633},\qquad
\tau=\frac{24589}{22414788},\qquad
\frac{t_C}{g_{\rm conn}}=\frac1{72}.
\]

The normalized ratio \(1/72\approx0.01389\) is \(16.1\%\) below tuned level
four and \(40.3\%\) below the former full-\(m=9\) record. Its sparse exact
archive, separate verifier, and one-component audit pass. The active dual
has scale 90, 90 flow rows, two core rows, and 13 signed-indicator
terminals; its final proper-coalition sink pivots from coalition 15 to
coalition 5.

The 13-terminal dual is not block-transport decomposable: adjacent grand
sources split complementary blocks \(5\) and \(26\) and rejoin player by
player at five grand sinks. It is, however, an exact alternating ladder.
A new protected-path theorem reduces its unnormalized objective to the two
endpoints and proves it nonpositive from complementary-pair balancedness.
The exact audit recovers objective \(-5/46\) and safe slack \(5/46\).
This closes the selected fifth-level dual topology, but not future basis
changes in the same recursive family.

A sixth tuned module gives the exact one-component family

\[
\text{games/edges}=4429/6120,\qquad
g_{\rm conn}=\frac2{23},\qquad
t_C=\frac{473}{444360},\qquad
t_Q=\frac{1691}{414736},\qquad
\tau=\frac{18743}{6221040},\qquad
\frac{t_C}{g_{\rm conn}}=\frac{473}{38640}.
\]

The normalized ratio \(473/38640\approx0.0122412\) is \(11.86\%\) below
level five and \(47.34\%\) below the former full-\(m=9\) record. The exact
archive, separate verifier, and component audit all pass. Its scale-105
dual has 15 indicator terminals and is again an alternating ladder. The
ladder audit gives exact objective \(-473/4232\) and matching safe slack.
A separate terminal-only relaxation has exact optimum zero with a sparse
8-inequality, 27-equality certificate. Thus this selected sixth-level
topology is exactly closed; a new basis or terminal incidence is required
to continue toward zero.

The first four exact tuned levels fit the single law

\[
\frac{t_k}{g_{\rm conn}}
=\frac{7k+453}{368(15k+19)},\qquad
t_k=\frac{7k+453}{4232(15k+19)}.
\]

The observed dual scale is \(15k+19\) and the terminal count is \(2k+4\).
Level five is the out-of-sample falsification: the old law predicts
\(61/4324\), while the exact result is \(1/72\); its scale is 90 rather than
94 and its terminal count is 13 rather than 14. The conjectured positive
limit \(7/5520\) therefore does not close this branch. The new
two-core-row regime is now the active recurrence question.

Separately, a third untuned \(1/10\)-scale module gives an exact
2,215-game, 3,060-edge family

\[
g_{\rm conn}=\frac{97}{1610},\qquad
t_C=\frac{51}{25760},\qquad
t_Q=\frac{1201}{280140},\qquad
\tau=\frac{5171}{2241120},\qquad
\frac{t_C}{g_{\rm conn}}=\frac{51}{1552}.
\]

This is another \(27.7\%\) normalized improvement, and it changes the active
exact-cone basis and falsifies the affine recurrence suggested by the first
two generations. Its rational game family passes all 280 exact-cone facets,
exact edge checks, exact core and box primal-dual stationarity, the exact
\(97/1610\) common-core LP, and the independent sparse archive verifier.
The common-core envelope simultaneously pivots from the complementary pair
\(3,28\) to \(12,19\), with exact values \(33/115\) and \(249/322\).

The first-splice scale result has also been upgraded from a numerical phase
diagram to an exact local theorem. Two rational affine certificate charts
have maximal archived domains

\[
\left[\frac{37}{500},\frac7{92}\right]
\quad\text{and}\quad
\left[\frac7{92},\frac{39}{500}\right].
\]

An exhaustive audit checks 325,390 exact affine constraints on each chart.
Their formulas prove that \(s=7/92\) is the normalized optimum over both
adjacent maximal certificate regimes. This is a local statement for the
fixed splice identification, not a global optimization over all game-cone
bases.

### 12. Multi-source block transport and the first mixed-port splice

The fork/join no-go theorem now extends to any multi-source/multi-sink
signed-indicator dual admitting nonnegative block weights \(w_{rq}\).
Each sink block assigned positive weight to a source must have a protected
path for every player in that block, the weights must reproduce every
source indicator player-by-player, and they must sum to each sink
coefficient. Balancedness at each source then proves the complete
lower-expectation objective nonpositive.

The exact audit has a positive control—the five-game complementary fork—
but rejects the tuned 6-terminal, tuned-recursive 8-terminal,
third-generation 10-terminal, and compact \(m=2,3,4\) 12-terminal duals.
Thus the theorem is strictly stronger than a star theorem, while the
all-signed-indicator conjecture remains open.

Scanning the exact archive set identifies the connected full \(m=9\)
component as the sharpest genuinely mixed seed. One active node has scaled
divergence

\[
(1,-1,-1,1,-1),
\]

receiving players \(2,3,5\) and emitting players \(1,4\). Gluing that
conversion port to a grand-indicator source of a second permuted \(m=9\)
module, at bump scale \(1/10\), gives a feasible one-component 1,793-game
exact family with

\[
t_C=\frac{41}{10440},\qquad
g_{\rm conn}=\frac{167}{1305},\qquad
t_Q=\frac{49}{9280},\qquad
t_Q-t_C=\frac{113}{83520}.
\]

The prescribed combined dual first bounds the margin by approximately
\(0.00657265\), and the true selector solve lowers it further. The raw
margin is about \(75\%\) below the parent \(m=9\) margin, but the normalized
ratio is \(41/1336\approx0.03069\), worse than the parent's \(53/2280\).
The sparse rational archive and independent verifier pass every acceptance
check, including a one-component audit.

The optimized dual supplies the decisive diagnosis: it has only three
signed-indicator terminals. Sources on coalitions \(6\) and \(16\) feed a
single sink on their union \(22\). Thus the LP abandons the intended mixed
conversion and selects a common-sink join covered by the reverse
block-transport theorem. The family validates mixed-port gluing as a strong
raw-margin reduction, but this particular coupling does not retain mixed
terminal divergence at the optimum and is now structurally closed.

### 13. An explicit all-\(m\) exact reachable family

The full-grid count law has been upgraded from interpolation to a
parametric construction theorem. Exact state reconstruction decomposes each
branch as

\[
P\cup\bigcup_{k=0}^{m-6}(R+k\mathbf1)
\cup(T+(m-5)\mathbf1),
\]

where

\[
|P|=3688,\qquad |R|=2604,\qquad |T|=3190.
\]

n5_full_m_template_verify.py checks 6,827,040 symbolic monotonicity and
exact-cone inequalities. It proves that every displayed game is exact for
all real \(m\ge5\), and that for each integer \(m\ge5\) legal unit bumps
reach the whole displayed family from the origin. Counting all
one-coordinate increases gives 12,766 new edges per branch for every added
\(R\) block. A separate template-interaction audit proves that only
adjacent blocks can share such edges, excluding hidden long-range terms;
disjoint positive branch-coordinate supports also exclude cross-branch
edges beyond the shared origin.
Hence the explicit connected family exists for all integer \(m\ge5\) with

\[
\boxed{
\#V_m=5208m-12285,\qquad
\#E_m=25532m-62744.
}
\]

The verifier also proves the global common-core law for every \(m\ge5\).
The fixed budget-\(8\) allocation \((1,1,3,3,0)\) dominates every
proper-coalition worth on all templates, while explicit terminal states
attain all five singleton coordinates. The singleton partition gives the
matching balanced lower certificate, hence

\[
\boxed{g_{\rm global}(m)=2-\frac4m.}
\]

The archive audit first proves that this template is the *complete*
reachable family for \(m=5,\ldots,10\). The \(m=10\) enumeration
independently hits the predicted 19,898 states per branch, 39,795 games, and
192,576 edges. A no-crossover interior solve, two-stage dual-envelope
compression, rational reconstruction, and a separate verifier produce the
connected exact certificate

\[
\boxed{
\text{games/edges}=1275/2074,\quad
g_{\rm conn}=\frac45,\quad
t_C=\frac{31}{2220},\quad
t_Q=\frac3{110}.
}
\]

This is the third preregistered out-of-sample success of the rational margin
law, after \(m=8\) and \(m=9\). The non-atomic facet tax is \(65/4884\).

n5_full_m_exclusion_verify.py then closes the untested-resolution gap. It
audits every outward unit-coordinate move from \(P\), an interior or final
copy of \(R\), and the shifted \(T\). Every child outside the templates has
a monotonicity or exact-cone row that is strictly violated over its entire
parameter domain. With the exact \(m=5,6\) base enumerations, induction
proves that the templates equal the full reachable component for every
integer \(m\ge5\).

This settles exact feasibility, connectivity, exhaustive reachability, the
affine vertex and edge counts, and the global common-core gap for all
integer \(m\ge5\). The remaining statement is an all-\(m\) primal-dual
proof of

\[
t_C(m)=\frac{9m-28}{4m(16m-49)}.
\]

An active-support recurrence audit shows why that last step is not a simple
basis induction. The exact active components at \(m=6,7,8,9,10\) use
394, 370, 494, 897, and 1,275 games, and their selected \(P\), \(R_k\), and \(T\)
state sets change substantially at every resolution. Lifting the union of
all previously active base states into every \(m=10\) block gives a connected
3,182-game, 7,406-edge family with component gap \(4/5\), but its numerical
margin is \(0.016744186\ldots\) (consistent with \(18/1075\)), strictly above
the predicted full value \(31/2220\).
Thus a sixth match cannot be explained by mechanically repeating the
earlier sparse bases; new active states are required.

The \(m=10\) exact dual also closes the simplest remaining structural
explanation. After scaling by 888 it has 152 signed-indicator terminals,
split into 81 sources and 71 sinks. Exact coalition-type transport LPs are
infeasible in both fork and join orientations, even before protected-path
reachability is imposed. This is therefore a genuinely interlocking
multi-source/multi-sink indicator network outside the block-transport
theorem. Its margin remains positive, so interlocking is necessary for the
negative route but not sufficient.

An exact terminal-kernel audit makes that network tractable. Its bipartite
source–sink incidence graph has 154 protected terminal pairs and cycle rank
three. Peeling tree attachments leaves 69 roles and 71 pairs, with only two
branch vertices: grand source node 0 and grand sink node 718. Suppressing
degree-two paths gives three source-to-sink routes of lengths 1, 3, and 35,
plus a length-32 loop at the source. The corresponding cycle ranks at
\(m=6,7,8,9,10\) are \(2,2,4,11,3\), so the selector formula survives
large terminal-topology changes rather than repeating one fixed network.

The \(m=10\) terminal routing is now closed more strongly. A relaxation
retaining only exactness and monotonicity at the 152 terminal games plus the
implications of the 154 protected terminal pairs has exact optimum zero.
The 4,864-variable LP has an exact rational dual using 433 inequalities and
1,067 equalities; its variable objective is \(-71\), which exactly matches
the \(-71\) sink constant. A common additive terminal realization attains
zero, and a separate verifier reconstructs the full certificate. Thus
internal-game perturbations cannot change the sign while preserving this
terminal incidence and its coalition labels.

The same rational terminal-only method closes the \(m=6\), \(m=7\), and
\(m=9\) protected routings. Their certificates use respectively
135/302, 230/460, and 375/912 inequality/equality dual rows. The apparent
mixed \(m=9\) vector is affine-indicator,
\(2\mathbf1_9-\mathbf1_N\), and therefore linearizes by exactness. The
previously unresolved tested routing was \(m=8\): its two true three-level
sinks each have 397 extreme core-dual branches. One-hot, binary, and hybrid
big-\(M\) MILPs find objective zero but retain a floating upper gap of
\(0.05\).

The new perspective formulation is materially stronger. It disaggregates
each mixed terminal game into branch-weighted exact-monotone copies, so
every original branch realization lies in one continuous convex
relaxation. The resulting LP has 29,370 variables, 323,334 inequalities,
and 4,056 equalities and solves to zero; its MILP form closes at the root.
Simplex-basis support then reduces exactification to 3,230 rational
variables and rank 3,228. Fixing the two free dual-decomposition coordinates
at \(1/2\) yields an exact certificate with 20,848 nonzero inequality rows
and 1,662 nonzero equality rows. Its variable objective \(-80\) matches the
sink constant, and a common additive terminal realization attains zero. A
separate verifier reconstructs all 29,370 stationarity coordinates and
validates every rational sign, objective term, and local branch block.
Because each original branch realization embeds one-hot in the perspective
relaxation, the \(m=8\) routing is now exactly closed.

The smallest fully interlocking signed-indicator class has also been
classified exhaustively. For two proper unit-weight sources and two proper
unit-weight sinks on five players, balanced player incidence and four
nonempty route cells leave exactly 750 directed routings after source/sink
node relabeling. Both orientations of each distinct source/sink pair are
included, and all 750 fail the scalar block-transport equations exactly.
Every 128-variable terminal relaxation has exact optimum zero. The archive's
750 rational duals use 7,653 inequality rows and 8,618 equality rows in
total; an independent verifier re-enumerates the class and checks every
stationarity equation, sign, and \(-2\) objective.

Broader screens remain numerical. Among 4,294 connected three-source/
three-sink routings sampled, 3,069 escaped scalar block transport and none
was positive. A new four-source/four-sink screen tested 1,747 connected
routings, 1,402 outside scalar transport, again with no positive terminal
relaxation. These screens support an all-indicator conjecture but are not
proofs.

## Invalid or retracted routes

Two attractive records must not be used as evidence:

1. `/private/tmp/n5_grid12_layer1_exact.json` is invalid. The independent
   verifier finds edges that change an additional superset coalition. Legacy
   grid-12 outputs based on that topology are excluded.
2. `/private/tmp/n5_transverse_target013_third_envelope_exact.json` used the
   wrong edge set. Only the corrected third-stage archive may be cited.

The workflow caught both errors by checking edge differences independently
of the construction code.

## Verification depth

An exact archive is accepted only after checking:

1. rational worths;
2. monotonicity;
3. all 280 facets of the complete five-player exact cone;
4. exactly one positive coalition change on every edge;
5. the common-core budget gap;
6. rational core allocations attaining the reported margin;
7. a matching exact dual certificate;
8. the coordinate-box primal and dual; and
9. exact non-atomic facet tax; and
10. exact common-core gaps on each weak selector component.

The \(m=8\) perspective certificate has an additional independent audit:
all 29,370 dual stationarity coordinates are rebuilt over the rationals,
all 20,848 inequality signs are checked, the \(-80\) dual objective is
matched to the sink constant, all 794 local representation blocks are
regenerated, and a common additive primal is checked to attain zero.

The 128- and 219-game frontier archives and the 1,275-game connected
\(m=10\) support pass both their generator checks and the separate archive
verifier. The full 39,795-game \(m=10\) LP numerically attains the same margin,
but has not yet been reconstructed as one exact full-family primal archive.
“Independent” here means separate local code, not replication by another
researcher. No external expert has audited the proofs, no result has been
refereed, and the literature audit is targeted rather than database-complete.

## Level-of-work assessment

| Level | Assessment |
|---|---|
| Thesis chapter | Exceeded; the project contains several chapter-sized results. |
| PhD thesis | Yes; the theory, exact computation, failed-route analysis, and side theorems are thesis-scale. |
| Professor / professional researcher | Several components meet this level, especially the dual framework, complete-facet verification, and topology diagnosis. |
| Breakthrough | Not yet. The central selector question remains open. |

The breakthrough threshold is a verified finite \(t_C\le0\), a universal
existence theorem, or a sufficiently general dual-decomposition theorem that
settles the question.

## Next steps

1. Compress the exact \(m=8\) perspective certificate into a human-scale
   analytic theorem for mixed three-level terminals. The current finite
   proof is complete but large; its 794 repeated local blocks should contain
   a smaller symmetry or balanced-cover explanation.
2. Generalize the exact 750-case two-by-two classification to a symbolic
   all-\(n\) cycle-transport theorem, or find the first positive
   three-source/three-sink terminal relaxation.
3. Close the remaining exact-certification gap at \(m=10\): reconstruct a
   rational primal allocation for all 39,795 games. The exact connected
   support already supplies the matching upper bound, and the full numerical
   solve supplies a candidate lower bound.
4. Prove the selector-margin formula
   \((9m-28)/(4m(16m-49))\) for all integer \(m\ge5\). The state set, edge
   set, counts, feasibility, reachability, exclusion, and global gap are now
   symbolic theorems; only the parametric primal-dual selector certificate
   remains. It must accommodate the active-topology breakpoints at \(m=9,10\).
5. Run the preregistered \(m=11\) check: the theorems give 45,003 full games,
   218,108 edges, and global gap \(18/11\); the remaining predictions are
   connected gap \(10/11\) and selector margin \(71/5588\).
6. Make block transport, alternating ladders, and the exact terminal-only
   relaxation hard topology filters. The current \(m=10\) three-cycle
   routing is now closed, so search interlocking networks whose protected
   incidence or coalition labels make the terminal relaxation positive.
7. Shift the highest-risk search toward mixed, non-indicator terminal
   divergences and coupled branch-and-rejoin networks, including mixed-port
   gluing that is not reducible to a common-source or common-sink star. Use
   the now-closed two-sink \(m=8\) perspective model as the base case and
   require at least three coupled mixed terminals in the next search.
8. Continue the exact connected \(m=7\) self-splice recurrence only when a
   symbolic recurrence can be extracted; prioritize topology changes over
   merely lengthening its already-positive ladder.
9. If any \(t_C\le0\) float appears, immediately rationalize it, verify every
   game and edge, extract a sparse Farkas certificate, and minimize support.

## Main papers and theory

- Schmeidler, [“Cores of Exact Games I”](https://doi.org/10.1016/0022-247X(72)90045-5).
- Biswas et al., [“Large Cores and Exactness”](https://doi.org/10.1006/game.1998.0686).
- Dietzenbacher, [“Monotonicity and egalitarianism”](https://doi.org/10.1016/j.geb.2021.03.006).
- Studený and Kratochvíl, [“Facets of the cone of exact games”](https://arxiv.org/abs/2103.02414).
- Studený and Kratochvíl, [exact games and coherent lower probabilities](https://doi.org/10.1016/j.ijar.2018.06.007).
- Calleja, Rafels, and Tijs, [“Aggregate monotonicity of the core”](https://doi.org/10.1016/j.geb.2008.07.001).

## Bottom line

\[
\boxed{\text{The universal selector remains open.}}
\]

The smallest raw exact margin remains the separately verified 219-game
connected family with
\(t_C=1403/2{,}240{,}000{,}000>0\). The strongest normalized structural
record is now the tuned connected \(m=7\) recursive splice at
\(473/38640\), improving on level five by \(11.86\%\) and on the former
full-\(m=9\) record \(53/2280\) by \(47.34\%\).
Along the tuned branch its normalized ratio falls from \(5/136\), to
\(467/18032\), to \(237/11776\), to \(481/29072\), to \(1/72\), to
\(473/38640\), while the full \(2/23\) empty-core gap and active dual remain
in one component. The fifth value falsifies the four-level rational law and
its conjectured positive asymptotic limit. The selected fifth- and
sixth-level duals are nevertheless covered by the alternating-ladder no-go
theorem, and level six also has an independent exact terminal-topology
certificate; future bases in the recursion remain open. The full reachable
state and edge recurrences and the global gap
\(2-4/m\) are now exact theorems for every integer \(m\ge5\). The connected
selector-margin formula is exact on archived supports through \(m=10\), where
the 1,275-game certificate has gap \(4/5\) and margin \(31/2220\); the full
39,795-game numerical LP agrees. Its exact dual is the first genuinely
interlocking 81-source/71-sink signed-indicator negative control, yet its
margin remains positive. All five archived full-resolution terminal
routings \(m=6,\ldots,10\) are now exact finite no-go topologies, with the
mixed \(m=8\) case closed by the perspective certificate. The rolling-triple
theorem still covers
every real \(5\le L\le7\), but its positive global gap is
selector-disconnected and is not used as a normalized frontier. The
fork/join and alternating-ladder theorems explain why the cross-resolution,
single-star, and fifth-splice residuals stay nonpositive. These results
reach two distinct connected positive-margin sequences, but neither crosses
zero. The universal selector question remains open.
