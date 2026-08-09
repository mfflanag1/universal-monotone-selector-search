# Alternating protected-cycle no-go theorem

**Status:** exact theorem  
**Scope:** arbitrary player count, arbitrary positive edge weights, arbitrary grand worths  
**Central selector question:** still open

## Result

Every alternating protected block-flow cycle has nonpositive total core
lower expectation.  The blocks may overlap arbitrarily, the edge weights
need not agree, and grand worths may increase along the protected paths.

Consequently, after protected player paths have been grouped into weighted
block arcs, every pseudoforest is harmless: each connected component may
have at most one independent undirected cycle.  Repeated leaf contraction
removes every tree attached to a unicyclic core, and the remaining alternating
cycle is covered here.

This strictly extends the previous computational six-cycle evidence.  It is
dimension-free and does not depend on the five- or six-player facet
catalogues.

## Setup

Index sources and sinks cyclically by

\[
r_0,\ldots,r_{m-1},\qquad q_0,\ldots,q_{m-1}.
\]

Source \(r_j\) has two outgoing protected block paths:

\[
r_j\longrightarrow q_j
  \quad\hbox{carrying }a_j\mathbf 1_{A_j},
\]

and

\[
r_j\longrightarrow q_{j-1}
  \quad\hbox{carrying }b_j\mathbf 1_{B_j},
\]

where indices are modulo \(m\) and all \(a_j,b_j>0\).  Thus

\[
d_{r_j}=a_j\mathbf1_{A_j}+b_j\mathbf1_{B_j},
\]

while sink \(q_j\) receives the matched \(A_j\)-arc and the preceding
\(B_{j+1}\)-arc:

\[
d_{q_j}=-a_j\mathbf1_{A_j}-b_{j+1}\mathbf1_{B_{j+1}}.
\]

Write the grand worths as \(G_{r_j}\) and \(G_{q_j}\).  Directed game
monotonicity gives \(G_{r_j}\le G_{q_j}\) along the matched path.  Write

\[
\underline E_v(z)=\min_{x\in C(v)}z\cdot x.
\]

## Two exactness bounds

For every exact game \(v\), coalitions \(A,B\), and \(a,b\ge0\),

\[
\boxed{
\underline E_v(a\mathbf1_A-b\mathbf1_B)
\le
(a-b)G+bv(N\setminus B)-av(N\setminus A).
}
\]

To prove it, exactness supplies \(x\in C(v)\) with

\[
x(N\setminus B)=v(N\setminus B).
\]

Efficiency and the identity

\[
a x(A)-b x(B)
=(a-b)G+b x(N\setminus B)-a x(N\setminus A)
\]

then give the displayed upper bound because
\(x(N\setminus A)\ge v(N\setminus A)\).

The unequal coefficients are important: the lemma is not restricted to a
signed indicator or to unit flow.

The cycle proof uses the following direct paired form.  Exactness at \(B\)
gives

\[
\underline E_r(a\mathbf1_A+b\mathbf1_B)
\le a[G_r-v^r(N\setminus A)]+bv^r(B).
\]

Indeed, choose a source-core allocation tight at \(B\), and use the
complementary core constraint to upper-bound its payoff to \(A\).
Exactness at \(N\setminus A\) similarly gives

\[
\underline E_q(-a\mathbf1_A-c\mathbf1_C)
\le-a[G_q-v^q(N\setminus A)]-cv^q(C).
\]

Here a sink-core allocation tight at \(N\setminus A\) maximizes its payoff
to \(A\), while its payoff to \(C\) is at least \(v^q(C)\).

## Cycle proof

Apply the paired bounds at \(r_j,q_j\) with the matched block \(A_j\), the
unmatched outgoing block \(B_j\), and the unmatched incoming block
\(B_{j+1}\).  Protected-path invariance along the matched path gives

\[
v^{r_j}(N\setminus A_j)=v^{q_j}(N\setminus A_j).
\]

Therefore

\[
\begin{aligned}
\underline E_{r_j}(d_{r_j})+
\underline E_{q_j}(d_{q_j})
\le{}& b_jv^{r_j}(B_j)
-b_{j+1}v^{q_j}(B_{j+1})\\
&+a_j(G_{r_j}-G_{q_j}).
\end{aligned}
\]

The grand term is nonpositive.  The path carrying \(B_j\) from \(r_j\) to
\(q_{j-1}\) gives

\[
v^{r_j}(B_j)\le v^{q_{j-1}}(B_j).
\]

Hence every positive block-worth term is at most the matching negative term
at the preceding sink.  Summing over \(j\), all block-worth terms telescope,
while every grand-worth difference is nonpositive.

It follows that

\[
\boxed{
\sum_j\underline E_{r_j}(d_{r_j})
+\sum_j\underline E_{q_j}(d_{q_j})\le0.
}
\]

## Consequences for the search

This theorem rules out all of the following as counterexample mechanisms:

- the unweighted alternating six-cycle;
- the weighted alternating six-cycle;
- every longer alternating source/sink cycle;
- arbitrary overlaps among the protected blocks;
- arbitrary positive flow weights; and
- every protected block-flow pseudoforest, including components with trees
  attached to a cycle and non-product networks with grand-coalition changes.

The theorem also explains why numerical searches over random five- and
six-player cycles repeatedly returned zero: the zero was structural, not an
optimization failure.

The highest-value remaining computational and theoretical target is now a
protected block-flow component of cycle rank at least two.  After its leaves
are removed, such a component necessarily has terminal degree at least three.
The minimal complete-bipartite example is \(K_{2,3}\).  Unequal grand worths
alone are not an escape hatch.

This result does not prove the universal coalitionally monotone selector.
