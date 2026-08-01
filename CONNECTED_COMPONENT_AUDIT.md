# Selector-component common-core audit

**Date:** 2026-07-30  
**Central selector question:** open

## Why this audit is necessary

Coalitional-monotonicity constraints act only along legal bump edges. If the
game graph has several weakly connected components, the selector problem
decomposes across them. A common-core allocation inside each component gives
a compatible constant selection on that component, even when the
intersection of the cores of *all* games in all components is empty.

Therefore a positive global common-core gap is not, by itself, evidence that
the component controlling the selector margin has an empty common core. A
counterexample must have \(t_C\le0\) in some component; a meaningful
empty-core frontier should report the common-core gap of that same component.

`n5_component_gap_audit.py` computes exact common-core gaps separately on
every weak selector component and identifies the component carrying the
active dual flow.

## Results

### Connected frontier that survives

The 219-game two-sided complement-path archive is one connected component:

\[
\text{games/edges}=219/218,\qquad
g_{\rm conn}=\frac{37}{140{,}000{,}000},\qquad
t_C=\frac{1403}{2{,}240{,}000{,}000}>0.
\]

Its global and component gaps are identical.

The 128-game one-sided predecessor, 310- and 356-game corrected transverse
archives, and 316- and 322-game hybrid archives are also single connected
components. Their previously reported global gaps are therefore genuine
selector-component gaps.

The full reachable-grid supports have the following active selector
components:

| \(m\) | Connected games / edges | Component gap | \(t_C\) |
|---:|---:|---:|---:|
| 5 | 195 / 291 | \(-2/15\) | \(17/620\) |
| 6 | 394 / 584 | \(0\) | \(13/564\) |
| 7 | 370 / 510 | \(2/7\) | \(5/252\) |
| 8 | 494 / 679 | \(1/2\) | \(11/632\) |
| 9 | 897 / 1,409 | \(2/3\) | \(53/3420\) |
| 10 | 1,275 / 2,074 | \(4/5\) | \(31/2220\) |

The \(m=6,7,8,9,10\) component gaps satisfy

\[
\boxed{
g_{\rm conn}(m)=2-\frac{12}{m}.
}
\]

Thus \(m=6\) is an exact connected common-core boundary, and every tested
\(m>6\) is a genuine connected empty-common-core family. Combining this
with the independently verified margin law gives the connected normalized
candidate

\[
\frac{t_C(m)}{g_{\rm conn}(m)}
=
\frac{9m-28}{8(m-6)(16m-49)}
\sim\frac{9}{128m}>0.
\]

This is the correct normalized weak-boundary sequence. The formula is
independently exactified through \(m=10\), not yet proved for all \(m\).

The uncompressed common-grand-\(7\) cross-resolution families also survive.
Their active selector components are:

| Resolutions | Connected games / edges | Component gap | \(t_C\) |
|---|---:|---:|---:|
| \(2,\ldots,5\) | 297 / 496 | \(1\) | \(1/36\) |
| \(2,\ldots,6\) | 415 / 632 | \(1\) | \(1/39\) |
| \(2,\ldots,7\) | 553 / 786 | \(1\) | \(1/42\) |

Their small compressed and player-swapped descendants do not preserve this
component gap; the full connected components are the accepted certificates.

### Constructions downgraded by the audit

The following archives have positive *global* gaps but nonpositive gaps in
their active selector components:

| Construction | Global gap | Active-component gap |
|---|---:|---:|
| \(m=2,3,4\) compact backbone | \(1\) | \(-1/3\) |
| Compressed \(m=2,\ldots,7\) best swap | \(1\) | \(-1\) |
| Rolling \(m=5,6,7\) best swap | \(6/5\) | \(-4/5\) |
| Rolling \(m=6,7,8\) best swap | \(4/3\) | \(-2/3\) |
| Rolling \(m=7,8,9\) best swap | \(10/7\) | \(-4/7\) |
| Full \(m=5,6\) cross-resolution support | \(6/5\) | \(-4/5\) |
| Blockwise addition through \(m=7\) | \(6/5\) | \(-4/5\) |
| Blockwise addition through \(m=8\) | \(6/5\) | \(-4/5\) |

These are valid exact margin calculations and may still illuminate dual
topology. They are not connected empty-common-core obstructions. Each active
component admits a constant common-core selection.

The symbolic rolling-triple theorem remains mathematically correct, but its
empty global common core is supplied across disconnected components. It
must be treated as a parametric exact-LP side theorem rather than evidence
that a connected selector obstruction approaches zero.

### Recursive self-splice diagnosis

Recursive splices of the \(m=2,3,4\) support produced exact margins as small
as \(67/798000\) while retaining a global gap \(27/14\). Component auditing
shows that their active selector components have gaps

\[
-\frac{17}{420},\quad
-\frac{17}{700},\quad
-\frac{11}{200},\quad
-\frac1{280},
\]

at successive audited stages. The shrinking displayed ratios were caused by
rescaling the active component while isolated envelope games maintained the
global gap. The sparse dual also collapsed first to balanced and then to
complementary indicator forks covered by `COMPLEMENTARY_FLOW_NO_GO.md`.

This recursive branch is therefore closed as a counterexample route.

### Connected \(m=7\) self-splices survive the audit

The separate self-splices of the full connected \(m=7\) component do not
have this defect:

| Family | Games / edges | Components | Active-component gap | Normalized margin |
|---|---:|---:|---:|---:|
| Base \(m=7\) | \(370/510\) | \(1\) | \(2/7\) | \(5/72\) |
| First splice | \(739/1020\) | \(1\) | \(9/161\) | \(197/3060\) |
| Second splice | \(1477/2040\) | \(1\) | \(9/161\) | \(401/8820\) |
| Tuned first splice | \(739/1020\) | \(1\) | \(2/23\) | \(5/136\) |
| Tuned second splice | \(1477/2040\) | \(1\) | \(2/23\) | \(467/18032\) |
| Tuned third splice | \(2215/3060\) | \(1\) | \(2/23\) | \(237/11776\) |
| Tuned fourth splice | \(2953/4080\) | \(1\) | \(2/23\) | \(481/29072\) |
| Tuned fifth splice | \(3691/5100\) | \(1\) | \(2/23\) | \(1/72\) |
| Tuned sixth splice | \(4429/6120\) | \(1\) | \(2/23\) | \(473/38640\) |
| Third splice | \(2215/3060\) | \(1\) | \(97/1610\) | \(51/1552\) |

For every displayed splice archive, the component containing the active
dual is the whole graph. Independent component LPs recover the displayed
exact gaps. The recursive reductions are therefore genuine connected
effects, not normalization by isolated envelope games.

## Research consequence

All future searches should enforce or audit at least one of:

1. a connected selector graph with positive common-core gap;
2. a positive common-core gap in the component carrying the active dual; or
3. direct incompatibility \(t_C\le0\), which needs no empty-core proxy.

Global gap divided by global margin must not be used as a strength metric
when the envelope and active dual live in different components.
