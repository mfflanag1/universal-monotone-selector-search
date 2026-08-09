# Universal monotone selector search — July 29 continuation

This is the focused continuation workspace for the universal coalitionally
monotone core-selector problem on exact cooperative games.

**Status (2026-08-01): open.** Expanded searches over interlocking flows,
permutation-orbit couplings, the full signed exact-game tangent cone, and
multilevel lattices remained nonnegative. The best margins approach zero only
through a degenerate additive-game limit; no counterexample or universal proof
has been obtained. New exact results classify every five-player three-level
mixed divergence, rule out the first four-terminal cycle with four maximal
local premiums, close every alternating protected cycle in arbitrary
dimension with arbitrary positive edge weights and arbitrary grand worths,
and prove a common-core theorem for all 20 irreducible maximally sharp
five-player \(K_{2,3}\) orbits at common grand worth.

## Orientation

- [`CURRENT_RESEARCH_UPDATE_2026-08-01.md`](CURRENT_RESEARCH_UPDATE_2026-08-01.md)
  is the latest status and verification summary.
- [`CURRENT_RESEARCH_REVIEW_2026-07-30.md`](CURRENT_RESEARCH_REVIEW_2026-07-30.md)
  records the preceding full-session review.
- [`COMPLEMENTARY_FLOW_NO_GO.md`](COMPLEMENTARY_FLOW_NO_GO.md) records the main
  structural obstruction.
- [`N5_MIXED_DIVERGENCE_PREMIUM_RESULT.md`](N5_MIXED_DIVERGENCE_PREMIUM_RESULT.md)
  gives the exact sharp premium table and the zero-premium 2–1–2 identity.
- [`N5_JOHNSON_SQUARE_NO_GO_RESULT.md`](N5_JOHNSON_SQUARE_NO_GO_RESULT.md)
  proves the exact 157,609-branch Johnson-square/cross-cosingleton no-go.
- [`ALTERNATING_CYCLE_NO_GO_RESULT.md`](ALTERNATING_CYCLE_NO_GO_RESULT.md)
  proves the dimension-free weighted alternating-cycle theorem and closes
  every protected block-flow pseudoforest, including unequal grand worths.
- [`N5_K23_COMMON_CORE_RESULT.md`](N5_K23_COMMON_CORE_RESULT.md) proves that
  all 20 irreducible maximally sharp cycle-rank-two \(K_{2,3}\) symmetry
  orbits have a common core for every common-grand realization, using 66,108
  exact rational envelope certificates.
- [`CONNECTED_GRID7_SELF_SPLICE_RESULT.md`](CONNECTED_GRID7_SELF_SPLICE_RESULT.md)
  and [`N5_FULL_REACHABLE_BOUNDARY_RESULT.md`](N5_FULL_REACHABLE_BOUNDARY_RESULT.md)
  document the strongest exact families and boundary certificates.

The consolidated project, including broader selector counterexamples and the
main claim ledger, is
[`mfflanag1/exact-game-monotone-selection`](https://github.com/mfflanag1/exact-game-monotone-selection).

## Highest-value review

The main question is whether the protected-flow inequality suggested by the
finite exact certificates can be proved for arbitrary legal comparison
networks, or whether a genuinely interlocking topology escapes it. The first
unclosed block-flow components have cycle rank at least two. At common grand
worth, the entire irreducible maximally sharp five-player \(K_{2,3}\) class is
now closed; unequal-grand or non-sharp \(K_{2,3}\) cases and more general
rank-two components remain. Additional forest or single-cycle brute-force
search is redundant.

Treat numerical margins as search diagnostics. Any claimed counterexample
must be reconstructed over the rationals and checked for exact-game facets,
legal one-coordinate edges, and a negative compatibility certificate.
