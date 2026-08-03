# Protected-path no-go for forks, joins, and alternating ladders

This note isolates a general structural reason that indicator-flow
forks and joins—including the recurring complementary pair/triple
constructions—stay compatible. The original balanced-cover theorem is
strictly strengthened here: the coalition blocks and their coverage may be
arbitrary.

## Setup

Consider a finite family of exact games and write \(G_v=v(N)\). Every
directed edge \(e:k\to\ell\) raises exactly one coalition \(S_e\), possibly
\(N\), and a protected flow for player \(i\) may use the edge only when
\(i\in S_e\).

For a game \(v\) and a coalition \(S\), exactness gives

\[
\underline E_v(\mathbf 1_S)=v(S),
\qquad
\underline E_v(-\mathbf 1_S)=v(N\setminus S)-G_v.
\]

## Protected-path invariance

**Lemma.** If player \(i\) has a directed protected path from game \(v\) to
game \(w\), then

\[
v(T)=w(T)
\]

for every coalition \(T\) that excludes \(i\).

**Proof.** Every edge on the path changes a coalition containing \(i\).
Because an edge changes exactly one worth coordinate, it leaves every
coalition excluding \(i\) unchanged. Iterate along the path. \(\square\)

This simple observation is stronger than ordinary reachability. Two games
can lie in the same undirected component while no protected player flow can
carry the dual divergence between them.

## Arbitrary-fork theorem

Let \(B_1,\ldots,B_m\) be arbitrary nonempty coalitions with nonnegative
weights \(a_1,\ldots,a_m\). Suppose a source game \(r\) sends the players of
\(B_j\), with flow mass \(a_j\), along protected paths to sink \(q_j\). The
divergences are

\[
d_r=\sum_j a_j\mathbf1_{B_j},
\qquad
d_{q_j}=-a_j\mathbf1_{B_j}.
\]

Put \(A=\sum_j a_j\). Exactness at each sink and protected-path invariance
give

\[
\underline E_{q_j}(-a_j\mathbf1_{B_j})
=a_j\left[v^r(N\setminus B_j)-G_{q_j}\right].
\]

For any \(x\in C(r)\),

\[
\begin{aligned}
\underline E_r\left(\sum_j a_j\mathbf1_{B_j}\right)
&\le \sum_j a_jx(B_j)\\
&=AG_r-\sum_j a_jx(N\setminus B_j)\\
&\le AG_r-\sum_j a_jv^r(N\setminus B_j).
\end{aligned}
\]

Directed paths imply \(G_{q_j}\ge G_r\). Adding the sink terms proves

\[
\sum_k\underline E_k(d_k)
\le\sum_j a_j(G_r-G_{q_j})
\le0.
\]

Therefore

\[
\boxed{\sum_k\underline E_k(d_k)\le0.}
\]

Thus no single-source indicator fork can certify incompatibility. No
balanced-cover assumption is needed.

## Arbitrary-join theorem

Conversely, suppose source \(q_j\) sends the players of \(B_j\), with mass
\(a_j\), along protected paths to a common sink \(u\). The divergences are

\[
d_{q_j}=a_j\mathbf1_{B_j},
\qquad
d_u=-\sum_j a_j\mathbf1_{B_j}.
\]

Exactness at each source and directed worth monotonicity give

\[
\sum_j\underline E_{q_j}(a_j\mathbf1_{B_j})
=\sum_j a_jv^{q_j}(B_j)
\le\sum_j a_jv^u(B_j).
\]

Meanwhile,

\[
\begin{aligned}
\underline E_u\left(-\sum_j a_j\mathbf1_{B_j}\right)
&=-\max_{x\in C(u)}\sum_j a_jx(B_j)\\
&\le-\sum_j a_jv^u(B_j).
\end{aligned}
\]

Hence every common-sink indicator join also has nonpositive objective:

\[
\boxed{\sum_k\underline E_k(d_k)\le0.}
\]

## Protected block-flow forest theorem

There is a broader decomposition that permits arbitrarily many mixed
transshipment nodes.  Replace every protected player flow by a source-to-sink
path decomposition.  Group paths with the same endpoints, mass, and player
set into abstract arcs

\[
e:u\longrightarrow v,qquad a_e\mathbf1_{B_e},quad a_e>0,
\]

where every player in \(B_e\) has a directed protected path from game \(u\)
to game \(v\).  Regard parallel block arcs as distinct edges.  If the
underlying undirected multigraph of these abstract arcs is a forest, then

\[
\boxed{\sum_k\underline E_k(d_k)\le0.}
\]

**Proof.**  Contract leaves.  If a leaf \(u\) is the source of its unique
incident block arc \(e:u\to v\), then

\[
\underline E_u(a_e\mathbf1_{B_e})
=a_ev^u(B_e)\le a_ev^v(B_e)
=\underline E_v(a_e\mathbf1_{B_e}).
\]

Let \(d_v\) be the old divergence at \(v\).  Deleting \(u\) and \(e\)
changes the divergence at \(v\) to
\(d'_v=d_v+a_e\mathbf1_{B_e}\).  Superadditivity of lower expectation gives

\[
\underline E_v(d'_v)
\ge \underline E_v(d_v)
  +\underline E_v(a_e\mathbf1_{B_e}),
\]

so replacing the leaf and its neighbor by the contracted neighbor can only
increase the objective.

If the unique incident arc is instead \(e:v\to u\), exactness,
protected-path invariance, and grand-worth monotonicity give

\[
\begin{aligned}
\underline E_u(-a_e\mathbf1_{B_e})
&=a_e\left[v^u(N\setminus B_e)-G_u\right]\\
&\le a_e\left[v^v(N\setminus B_e)-G_v\right]\\
&=\underline E_v(-a_e\mathbf1_{B_e}).
\end{aligned}
\]

Now deleting the leaf changes the neighbor divergence to
\(d'_v=d_v-a_e\mathbf1_{B_e}\), and the same superadditivity argument applies.
Iterating leaves reduces every tree component to one node with zero
divergence and lower expectation zero.  Since every contraction upper-bounds
the previous objective, the original objective is nonpositive. \(\square\)

The multigraph qualification matters.  Two differently weighted or
differently supported blocks with the same endpoints are parallel edges and
already create cycle rank one unless they can be combined into a single
indicator block.  Thus the theorem does not silently split an aggregated
terminal lower expectation into separate terms; that split would have the
wrong inequality direction.

This theorem puts a sharp topological condition on any counterexample:
after every possible grouping of equal-mass protected player paths, its
terminal path multigraph must contain an undirected cycle.  The exhaustive
two-source/two-sink class below is exactly the first fully interlocking
four-cycle beyond this forest theorem.

## Alternating protected-cycle theorem

The forest condition can be relaxed for source-to-sink block-flow graphs of
maximum degree two.  Every alternating cycle is also harmless, for arbitrary
player count, arbitrary overlapping blocks, arbitrary positive edge weights,
and arbitrary grand-worth increases along the protected paths.

Index sources and sinks cyclically.  Source \(r_j\) sends
\(a_j\mathbf1_{A_j}\) to its matched sink \(q_j\), and sends
\(b_j\mathbf1_{B_j}\) to \(q_{j-1}\).  Pair the source and sink terms across
the matched \(A_j\)-arc.  Exactness at \(B_j\) and at
\(N\setminus A_j\), followed by protected-path invariance, gives

\[
\begin{aligned}
\underline E_{r_j}(d_{r_j})+
\underline E_{q_j}(d_{q_j})
\le{}&b_jv^{r_j}(B_j)
-b_{j+1}v^{q_j}(B_{j+1})\\
&+a_j(G_{r_j}-G_{q_j}).
\end{aligned}
\]

The grand term is nonpositive.  Monotonicity along the unmatched
\(B_j\)-path bounds \(v^{r_j}(B_j)\) by
\(v^{q_{j-1}}(B_j)\), so the block-worth terms telescope around the cycle.
Therefore

\[
\boxed{\sum_k\underline E_k(d_k)\le0.}
\]

Together with leaf contraction from the forest theorem, this closes every
protected block-flow pseudoforest: each component may contain at most one
independent undirected cycle, with arbitrary trees attached.  A complete
proof is recorded in `ALTERNATING_CYCLE_NO_GO_RESULT.md`.

## Sharp five-player K2,3 common-core classification

The first irreducible cycle-rank-two class is also exactly safe at common
grand worth. Consider unit protected flows on a \(K_{2,3}\) with two
degree-three sources and three degree-two sinks. Require every source to have
all four incidence levels, every pair within either source row to have sharp
five-player premium \(1/2\), every sink to have sharp premium \(1/2\), the
union of the six blocks to be \(N\), and their intersection to be empty.

There are 4,440 configurations after normalizing the sink labels against the
first source row, and exactly 20 orbits after player, source, and sink
symmetries. For every orbit and every common-grand monotone exact realization,

\[
\boxed{\bigcap_{k=0}^4 C(v^k)\ne\varnothing.}
\]

The proof passes to the pointwise worth envelope of the three sink games.
Protected invariance leaves only two through six sink coordinates that can
differ. The five-player balanced-weight polytope has exactly 1,291 extreme
points, enumerated with integer determinants and Cramer's rule. Expanding the
remaining envelope maxima gives 66,108 linear cases across the 20 orbits.
Every case has an exact rational dual bound equal to the grand worth, with
zero reconstruction failures. Thus every sink envelope is balanced, its core
is the common intersection of the three sink cores, and sink dominance places
that same point in both source cores.

A common core point makes every zero-sum terminal divergence safe, not merely
the unit divergence used to identify the topology. The full statement and
reproduction commands are in `N5_K23_COMMON_CORE_RESULT.md`.

## Single mixed-hub theorem

The fork and join proofs are the two degenerate cases of a stronger local
result.  Let distinct indicator sources \(r_j\) feed one hub \(h\), and let
the hub feed distinct indicator sinks \(q_\ell\).  Source \(r_j\) transports
the block \(B_j\) with mass \(a_j>0\) to \(h\); the hub transports the block
\(C_\ell\) with mass \(b_\ell>0\) to \(q_\ell\).  Every transported player
has a directed protected path between the stated endpoints.  The
divergences are

\[
d_{r_j}=a_j\mathbf1_{B_j},\qquad
d_h=\sum_\ell b_\ell\mathbf1_{C_\ell}
      -\sum_j a_j\mathbf1_{B_j},\qquad
d_{q_\ell}=-b_\ell\mathbf1_{C_\ell}.
\]

Then, even though \(d_h\) may have several levels and mixed signs,

\[
\boxed{\sum_k\underline E_k(d_k)\le0.}
\]

Indeed, directed worth monotonicity from \(r_j\) to \(h\) gives

\[
\underline E_{r_j}(a_j\mathbf1_{B_j})
=a_jv^{r_j}(B_j)\le a_jv^h(B_j).
\]

Protected-path invariance and grand-worth monotonicity from \(h\) to
\(q_\ell\) give

\[
\begin{aligned}
\underline E_{q_\ell}(-b_\ell\mathbf1_{C_\ell})
&=b_\ell\left[v^{q_\ell}(N\setminus C_\ell)-G_{q_\ell}\right]\\
&\le b_\ell\left[v^h(N\setminus C_\ell)-G_h\right].
\end{aligned}
\]

Choose \(x\in C(h)\) minimizing \(d_h\cdot x\).  Core feasibility implies

\[
v^h(B_j)\le x(B_j),\qquad
v^h(N\setminus C_\ell)-G_h\le-x(C_\ell).
\]

Consequently all non-hub terms sum to at most \(-d_h\cdot x\), while the
hub term is exactly \(d_h\cdot x\).  Their total is nonpositive.

This is also the star-shaped case of the protected block-flow forest theorem.
The direct proof remains useful because it shows explicitly how one hub core
allocation absorbs an arbitrary mixed divergence.  A large non-additivity
premium at one mixed node cannot by itself produce a counterexample.  The
requirement that each terminal divergence be one signed indicator is
substantive, because lower expectation is superadditive and a terminal
carrying several parallel blocks cannot generally be split without loss.

## Complete-bipartite sink-overlap reduction

There is a dimension-free reduction for every two-source/two-sink protected
flow, even when all four terminal divergences are mixed.  Let sources
(s_0,s_1) feed both sinks (t_2,t_3), and let the four divergences sum to
zero.  Pointwise game monotonicity gives (s_0\le t_2) and (s_1\le t_2),
so monotonicity and superadditivity of lower expectation imply

\[
\begin{aligned}
\sum_k\underline E_k(d_k)
&\le \underline E_{t_2}(d_{s_0})
 +\underline E_{t_2}(d_{s_1})
 +\underline E_{t_2}(d_{t_2})
 +\underline E_{t_3}(d_{t_3})\\
&\le \underline E_{t_2}(-d_{t_3})
 +\underline E_{t_3}(d_{t_3}).
\end{aligned}
\]

The last expression is nonpositive exactly when the scalar ranges overlap:

\[
\min_{x\in C(t_2)}(-d_{t_3})\cdot x
\le
\max_{y\in C(t_3)}(-d_{t_3})\cdot y.
\]

Thus a (K_{2,2}) counterexample requires two incomparable sink cores whose
relevant one-dimensional projections are strictly separated in the wrong
order.  The symmetric reduction through (t_3) gives a second necessary
separation condition involving (d_{t_2}).  This replaces a four-game mixed
optimization by two paired-core interval tests.  The Johnson theorem below
proves the required overlap for the first maximally dangerous five-player
instance.

## Exact five-player Johnson-square theorem

The first maximally dangerous cycle beyond the forest theorem is also safe.
On players (0,1,2,3,4), take the four protected blocks

\[
012,\quad013,\quad024,\quad034
\]

on a directed (K_{2,2}) from sources (s_0,s_1) to sinks (t_2,t_3), in
row-major order.  The resulting divergences are

\[
\begin{aligned}
d_{s_0}&=(2,2,1,1,0),&d_{s_1}&=(2,0,1,1,2),\\
d_{t_2}&=(-2,-1,-2,0,-1),&d_{t_3}&=(-2,-1,0,-2,-1).
\end{aligned}
\]

All four have sharp local five-player premium (1/2), yet

\[
\boxed{\sum_k\underline E_k(d_k)\le0.}
\]

The protected equalities force the games to agree except at four
co-singletons: one relaxed coordinate private to each source and one tightened
coordinate private to each sink.  In particular, (t_2,t_3) differ only at
(A=N\setminus\{2\}) and (B=N\setminus\{3\}), with opposite directions.
For (e=(2,1,0,2,1)=-d_{t_3}), an exact exhaustive certificate proves the
cross-cosingleton swap inequality

\[
\underline E_{t_2}(e)+\underline E_{t_3}(-e)\le0.
\]

Indeed, each lower expectation has 397 extreme core-dual representations.
All (397^2=157{,}609) affine branch pairs have exact rational dual upper
bounds zero over the complete paired monotone exact-game polytope.  The
verifier found zero failures; floating residuals are not used as proof.

Finally, game monotonicity and lower-expectation superadditivity give

\[
\begin{aligned}
\sum_k\underline E_k(d_k)
&\le \underline E_{t_2}(d_{s_0})
 +\underline E_{t_2}(d_{s_1})
 +\underline E_{t_2}(d_{t_2})
 +\underline E_{t_3}(d_{t_3})\\
&\le \underline E_{t_2}(e)+\underline E_{t_3}(-e)\le0.
\end{aligned}
\]

The full derivation and reproduction command are in
`N5_JOHNSON_SQUARE_NO_GO_RESULT.md`; the exact run summary is
`results/n5_cross_cosingleton_branch_sweep_exact.json`.

## Block-transport theorem for several sources and sinks

The fork and join arguments extend to a genuinely multi-source/multi-sink
class.  Let the positive and negative terminal divergences be

\[
d_r=a_r\mathbf1_{B_r},\qquad
d_q=-b_q\mathbf1_{C_q},
\]

where all coefficients are positive.  Suppose there are nonnegative block
weights \(w_{rq}\) satisfying

\[
\sum_{q:i\in C_q}w_{rq}=a_r\mathbf1_{\{i\in B_r\}}
\quad(r,i),
\qquad
\sum_rw_{rq}=b_q
\quad(q).
\]

Also require that whenever \(w_{rq}>0\), every player \(i\in C_q\) has a
directed protected path from \(r\) to \(q\).  Then

\[
\boxed{\sum_k\underline E_k(d_k)\le0.}
\]

**Proof.** Split the coefficient of each sink term according to
\(b_q=\sum_rw_{rq}\). Protected-path invariance and grand-worth
monotonicity give, whenever \(w_{rq}>0\),

\[
v^q(N\setminus C_q)=v^r(N\setminus C_q),
\qquad G_q\ge G_r.
\]

The terms assigned to one source \(r\) are therefore at most

\[
a_rv^r(B_r)+\sum_qw_{rq}
\left[v^r(N\setminus C_q)-G_r\right].
\]

Put \(A_r=\sum_qw_{rq}\).  At game \(r\), the weighted collection containing
\(B_r\) with weight \(a_r\) and every \(N\setminus C_q\) with weight
\(w_{rq}\) is balanced with coverage \(A_r\). Indeed, a player in \(B_r\)
receives \(a_r+A_r-a_r=A_r\), while a player outside \(B_r\) receives
\(A_r-0=A_r\). Balancedness of \(r\) bounds the displayed expression by
zero. Summing over sources proves the claim. \(\square\)

There is a reverse, sink-grouped block transport. Suppose instead that
nonnegative weights \(z_{rq}\) satisfy

\[
\sum_qz_{rq}=a_r
\quad(r),
\qquad
\sum_{r:i\in B_r}z_{rq}
=b_q\mathbf1_{\{i\in C_q\}}
\quad(q,i),
\]

and whenever \(z_{rq}>0\), every player \(i\in B_r\) has a protected path
from \(r\) to \(q\). Split each source coefficient according to
\(a_r=\sum_qz_{rq}\). Directed monotonicity gives

\[
v^r(B_r)\le v^q(B_r).
\]

For one sink \(q\), the collection containing every \(B_r\) with weight
\(z_{rq}\), together with \(N\setminus C_q\) with weight \(b_q\), is
balanced with coverage \(b_q\). Hence

\[
\sum_rz_{rq}v^q(B_r)
+b_q\left[v^q(N\setminus C_q)-G_q\right]\le0.
\]

Summing over sinks proves the same nonpositivity result.

Together, the source-grouped and sink-grouped theorems are strictly broader
than a star: several sources and sinks may be active and terminal
coefficients may be split. Their scope is nevertheless exact. In the first
orientation a sink's whole player block must be transportable from one
source as a block; in the reverse orientation a source's whole block must
be transportable to one sink as a block.

## Exhaustive five-player two-by-two interlocking classification

The smallest fully interlocking unit-indicator class can be closed even when
both scalar block-transport orientations fail. Take two proper source blocks
\(A_1,A_2\), two proper sink blocks \(C_1,C_2\), and require balanced player
incidence

\[
\mathbf1_{A_1}+\mathbf1_{A_2}
=\mathbf1_{C_1}+\mathbf1_{C_2}.
\]

Route every player occurrence from a source containing that player to a sink
containing it, and require all four source-sink route cells to be nonempty.
Players lying in both source and both sink blocks admit two crossing
orientations.

After quotienting only by relabeling the two sources and the two sinks, there
are exactly 750 directed five-player routings with distinct source and sink
pairs. Both directions of every pair are included. Exact rank checks show
that all 750 lie outside scalar block transport.

For every routing, the 128-variable terminal relaxation has exact optimum
zero. The archive contains 750 rational dual certificates, with 7,653
nonzero inequality rows and 8,618 nonzero equality rows in total. Every dual
has variable objective \(-2\), matching the two sink constants, and a common
additive terminal game attains zero. The independent verifier:

1. re-enumerates all 750 routings;
2. checks exact player-incidence balance and route supports;
3. proves the scalar transport equations inconsistent over the rationals;
4. reconstructs every exact-cone and protected-path dual row; and
5. verifies stationarity, signs, and objective for all 750 certificates.

This is a finite five-player classification, not yet the sought
all-signed-indicator theorem. It does show that the first fully interlocking
class beyond block transport is exactly safe.

## Alternating-ladder theorem

There is a second genuinely multi-source/multi-sink class not covered by
block transport. Let \(A,B\) be a nontrivial partition of \(N\), and suppose
all games have common grand worth \(G\). Take grand-indicator sources
\(r_0,\ldots,r_k\), grand-indicator sinks
\(q_0,\ldots,q_{k-1}\), a \(B\)-indicator sink \(\ell\), and an
\(A\)-indicator sink \(u\). Protected player paths form the alternating
ladder

\[
\begin{array}{lll}
r_0:A\to q_0,&r_0:B\to\ell,\\
r_j:A\to q_j,&r_j:B\to q_{j-1}&(0<j<k),\\
r_k:A\to u,&r_k:B\to q_{k-1}.
\end{array}
\]

Then its lower-expectation objective is nonpositive.

Indeed, the grand source/sink terms telescope, while exactness and
protected-path invariance at the endpoints give

\[
\mathcal O=v^\ell(A)+v^u(B)-G
=v^{r_0}(A)+v^{r_k}(B)-G.
\]

At every intermediate sink \(q_j\), the \(A\)-path from \(r_j\) fixes the
worth of \(B\):

\[
v^{q_j}(B)=v^{r_j}(B).
\]

The directed \(B\)-path from \(r_{j+1}\) to the same sink and coordinatewise
worth monotonicity give

\[
v^{r_{j+1}}(B)\le v^{q_j}(B)=v^{r_j}(B).
\]

Thus \(v^{r_k}(B)\le v^{r_0}(B)\). Balancedness of the complementary pair
\(A,B\) at \(r_0\) now yields

\[
\boxed{
\mathcal O
\le v^{r_0}(A)+v^{r_0}(B)-G
\le0.
}
\]

This proof depends on player-split rejoining at the grand sinks, so it is
not a corollary of block transport.

### Exact decomposition audit

`n5_indicator_transport_audit.py` extracts signed-indicator terminals from
an exact margin dual, computes player-protected reachability, and solves the
rational feasibility system above.  The tests distinguish the proved class
from the unresolved residual class:

| Exact archive | Terminals | Block transport |
|---|---:|---:|
| Five-game complementary boundary | 3 | yes |
| Tuned first \(m=7\) splice | 6 | no |
| Tuned second \(m=7\) splice | 8 | no |
| Third-generation \(m=7\) splice | 10 | no |
| Tuned fifth \(m=7\) splice | 13 | no; alternating ladder |
| Compact \(m=2,3,4\) overlap | 12 | no |
| Exact \(m=9\) mixed-port splice | 3 | yes, join orientation |

The positive test decomposes the grand source into the complementary sink
blocks \(6\) and \(25\), both with weight one.  The negative tests are not
counterexamples to selector existence: their exact selector margins remain
strictly positive, so their lower-expectation residuals have the safe sign.
For the tuned first and third-generation splices, decomposition already
fails in the player-by-player block equations before protected reachability
is imposed. They show that a theorem for *all* signed-indicator networks, if
true, needs a stronger argument than block transport.

The mixed-port construction is a useful negative control. It was designed
from a genuine \((1,-1,-1,1,-1)\) conversion node, but after optimizing the
full selector LP its exact dual collapses to sources on coalitions \(6\)
and \(16\) feeding one sink on their union \(22\). The reverse
block-transport theorem certifies this selected join topology as
nonpositive.

The tuned fifth splice supplies a different positive control. Its six
grand sources split blocks \(5\) and \(26\) across five grand sinks and two
proper endpoint sinks. It fails both block-transport orientations because
each grand sink mixes players from adjacent sources. The alternating-ladder
theorem nevertheless applies. `n5_alternating_ladder_audit.py` verifies the
exact protected paths, worth invariances, and slack decomposition:

\[
\mathcal O=-\frac5{46},\qquad
-\mathcal O=\frac5{46}.
\]

### Exact terminal-cone closure of five full-resolution routings

The full \(m=10\) active dual is harder: its 152 signed-indicator terminals
form a connected incidence graph of cycle rank three, and both
block-transport orientations are algebraically infeasible. It is not an
alternating ladder.

`n5_terminal_topology_relaxation.py` removes every internal game and keeps
only necessary terminal conditions: the complete exact cone, game
monotonicity, directed worth monotonicity, and protected-path invariance.
The resulting 4,864-variable relaxation has exact optimum zero. Its rational
dual uses 433 inequalities and 1,067 equalities and proves variable
objective \(-71\), exactly matching the sink constant. The independent
`n5_verify_terminal_topology_certificate.py` rebuilds the rows and verifies
stationarity and objective over the rationals. A common additive terminal
game attains zero.

This is a finite no-go theorem for the archived terminal incidence and
coalition labels: no internal realization of those protected routes can
have positive lower-expectation objective. It is not a decomposition
theorem for arbitrary interlocking indicator networks.

The same terminal-only construction and independent verification close the
full \(m=6\), \(m=7\), and \(m=9\) routings at exact objective zero. The
sole mixed-sign \(m=9\) vector is affine-indicator,
\(2\mathbf1_9-\mathbf1_N\), so its lower expectation is
\(2v(9)-v(N)\). The \(m=8\) dual is excluded from this statement: its two
three-level sink vectors require a maximum over 397 core-dual branches
each. The original one-hot, binary, and hybrid big-\(M\) branch MILPs all
find a zero incumbent but retain upper gap \(0.05\).

A stronger perspective formulation disaggregates each mixed terminal game
over all 397 extreme core-dual branches. Every original branch realization
is feasible in the resulting continuous convex relaxation. Its
29,370-variable LP has 323,334 inequalities and 4,056 equalities and solves
to objective zero; the corresponding perspective MILP closes at the root.

Simplex-basis support reduces exactification to a 3,230-variable rational
system of rank 3,228. Fixing its two harmless dual-decomposition degrees of
freedom at \(1/2\) gives an exact certificate with 20,848 nonzero
inequality rows and 1,662 nonzero equality rows. Its variable objective is
\(-80\), exactly matching the \(-80\) sink constant. A common additive
terminal realization attains zero. The independent verifier rebuilds all
29,370 stationarity coordinates, checks every rational sign and the
objective, and validates every one of the 794 extreme-representation
blocks. Therefore the archived \(m=8\) routing is also an exact finite
no-go topology.

The remaining sections give earlier balanced-cover derivations and useful
specializations on a common-grand section. They are now corollaries of these
analytic theorems.

## Balanced-fork theorem

Let \(B_1,\ldots,B_m\) be nonempty player blocks with nonnegative weights
\(a_1,\ldots,a_m\). Suppose they form a balanced collection of coverage
\(c\):

\[
\sum_{j:i\in B_j}a_j=c
\qquad(i\in N).
\]

Suppose a source game \(r\) sends the players of \(B_j\), with flow mass
\(a_j\), along protected paths to sink \(q_j\). The divergences are

\[
d_r=c\mathbf1_N,
\qquad
d_{q_j}=-a_j\mathbf1_{B_j}.
\]

Put \(A=\sum_j a_j\). Exactness and protected-path invariance give

\[
\begin{aligned}
\sum_k\underline E_k(d_k)
&=
cG+\sum_j a_j
\left[v^{q_j}(N\setminus B_j)-G\right]\\
&=
\sum_j a_j v^r(N\setminus B_j)-(A-c)G.
\end{aligned}
\]

The complemented collection
\((N\setminus B_j,a_j)\) is balanced with coverage \(A-c\), because every
player is excluded from blocks of total weight \(A-c\). Balancedness of the
source game therefore gives

\[
\sum_j a_j v^r(N\setminus B_j)\le(A-c)G.
\]

Hence

\[
\boxed{\sum_k\underline E_k(d_k)\le0.}
\]

No single-source balanced fork of indicator flows can be an incompatibility
certificate.

## Balanced-join theorem

The reverse orientation is also harmless. Suppose source \(q_j\) sends
players \(B_j\), with mass \(a_j\), along directed protected paths to one
sink \(u\). The divergences are

\[
d_{q_j}=a_j\mathbf1_{B_j},
\qquad
d_u=-c\mathbf1_N.
\]

Worths are coordinatewise nondecreasing along every directed path, so

\[
v^{q_j}(B_j)\le v^u(B_j).
\]

Consequently,

\[
\begin{aligned}
\sum_k\underline E_k(d_k)
&=
\sum_j a_j v^{q_j}(B_j)-cG\\
&\le
\sum_j a_j v^u(B_j)-cG\\
&\le0,
\end{aligned}
\]

where the final inequality is the balancedness inequality at the common
sink.

Thus no common-sink balanced join of indicator flows can be an
incompatibility certificate.

## Subset-balanced extension

The source blocks need not cover all players. Let \(B\subseteq N\), and let
\((B_j,a_j)\) be a weighted cover of \(B\) with constant coverage \(c\):

\[
B_j\subseteq B,
\qquad
\sum_{j:i\in B_j}a_j=c
\quad(i\in B).
\]

### Join on a subset

Suppose source \(q_j\) sends block \(B_j\) to one sink \(u\). The sink
divergence is \(-c\mathbf1_B\). Exactness gives

\[
\underline E_u(-c\mathbf1_B)
=c\left[v^u(N\setminus B)-G\right].
\]

Directed monotonicity and balancedness at \(u\) imply

\[
\begin{aligned}
\sum_k\underline E_k(d_k)
&=
\sum_j a_jv^{q_j}(B_j)
+c\left[v^u(N\setminus B)-G\right]\\
&\le
\sum_j a_jv^u(B_j)
+c\,v^u(N\setminus B)-cG\\
&\le0.
\end{aligned}
\]

The final line uses the balanced collection consisting of the \(B_j\), with
their weights \(a_j\), plus \(N\setminus B\) with weight \(c\).

### Fork on a subset

Conversely, suppose source \(r\) sends block \(B_j\) to sink \(q_j\). Its
divergence is \(c\mathbf1_B\). Protected-path invariance gives

\[
v^{q_j}(N\setminus B_j)=v^r(N\setminus B_j).
\]

Writing \(A=\sum_j a_j\), the objective becomes

\[
c\,v^r(B)
+\sum_j a_jv^r(N\setminus B_j)-AG.
\]

The collection containing \(B\) with weight \(c\) and every
\(N\setminus B_j\) with weight \(a_j\) has coverage \(A\) on every player:
inside \(B\), the complemented blocks contribute \(A-c\), while outside
\(B\) they contribute \(A\). Balancedness again makes the objective
nonpositive.

The 322-game empty/taxed recursive hybrid is an exact example of the join
case. Its non-grand ports are sources on coalitions \(2\) and \(24\) and a
sink on their union \(26\); exactness at the sink contributes coalition
\(5=N\setminus26\). The source blocks and coalition 5 partition the player
set, so this dual mechanism cannot cross zero.

## Complementary fork as a corollary

Fix a nonempty proper coalition \(S\). Suppose a source game \(r\) sends:

- every player in \(N\setminus S\) along protected paths to sink \(p\); and
- every player in \(S\) along protected paths to sink \(q\).

Give every routed player flow mass \(a>0\). The three nonzero divergences are

\[
d_r=a\mathbf 1_N,\qquad
d_p=-a\mathbf 1_{N\setminus S},\qquad
d_q=-a\mathbf 1_S.
\]

Protected-path invariance gives

\[
v^p(S)=v^r(S),\qquad
v^q(N\setminus S)=v^r(N\setminus S).
\]

This is the balanced-fork theorem for the two blocks
\(B_1=N\setminus S\) and \(B_2=S\), both of weight \(a\). Directly, the
lower-expectation objective is

\[
\begin{aligned}
\underline E_r(d_r)
+\underline E_p(d_p)
+\underline E_q(d_q)
&=
a\left[
G+v^p(S)-G+v^q(N\setminus S)-G
\right]\\
&=
a\left[
v^r(S)+v^r(N\setminus S)-G
\right]\\
&\le0.
\end{aligned}
\]

The last inequality follows from core nonemptiness, since the two
complementary coalition constraints must fit inside the common budget.

Thus a complementary empty-envelope gap attained at two incomparable games
cannot be converted into a positive fork-flow objective: the protected paths
force the relevant worths back to one source game, where balancedness blocks
the violation.

## Complementary join as a corollary

The reverse orientation is also harmless. Suppose players in \(S\) flow from
source \(p\) to a common sink \(u\), while players in \(N\setminus S\) flow
from source \(q\) to \(u\). The objective is

\[
a\left[v^p(S)+v^q(N\setminus S)-G\right].
\]

Worths are coordinatewise nondecreasing along every directed path, so

\[
v^p(S)\le v^u(S),
\qquad
v^q(N\setminus S)\le v^u(N\setminus S).
\]

Core nonemptiness of \(v^u\) then gives

\[
v^p(S)+v^q(N\setminus S)
\le
v^u(S)+v^u(N\setminus S)
\le G.
\]

Hence the join-flow objective is nonpositive.

## Exact interlocking negative control at \(m=10\)

The exact connected \(m=10\) full-grid support shows that the theorem's
hypothesis is genuinely restrictive. Its core dual has 152 signed-indicator
terminals—81 sources and 71 sinks—all with unit coefficient after scaling by
888. Nevertheless it admits neither orientation of block transport.

`n5_indicator_transport_audit.py` first aggregates terminals by coalition
type. This quotient is exact: any node-level transport can be averaged over
equal source and sink types, and any quotient transport can be split back
proportionally among equal types. Exact rational feasibility LPs on the
quotient are infeasible for both the source-grouped fork equations and the
sink-grouped join equations. Reachability restrictions therefore cannot
restore either decomposition.

Thus the 1,275-game, 2,074-edge connected certificate with

\[
g_{\rm conn}=\frac45,qquad t_C=\frac{31}{2220}>0
\]

is a fully interlocking multi-source/multi-sink indicator network outside
the no-go theorem. It is still compatible, so escaping block transport is
necessary for this route but not sufficient for a counterexample.

## Research consequence

The result covers every single-source or single-sink indicator fork/join,
with arbitrary overlapping blocks and weights, and every multi-terminal
network admitting the block-transport decomposition above. This includes the
complementary mechanisms that recur in the bridge chains, player partitions,
the 322-game hybrid, co-singleton stars, and the new sharp-cap
cross-resolution dual. It explains why:

- moving complementary envelope maxima to different games is insufficient;
- serial bridges can enlarge a dual normalization denominator without
  changing the positive residual; and
- a direct common source or sink does not turn the empty common core into a
  monotonicity obstruction.

It does **not** prove the universal selector theorem or even the
all-signed-indicator conjecture. A remaining counterexample must avoid
decomposition into these pieces.  In particular, its reduced protected
block-flow component must have cycle rank at least two.  It could use:

- several interlocking sources and sinks with no star decomposition;
- mixed, non-indicator divergence vectors at multiple internal nodes;
- coupled branch-and-rejoin networks not decomposable into the fork/join or
  pseudoforest patterns above; or
- varying grand worths combined with one of those higher-rank topologies.

Those are the topology classes that future searches should target.
