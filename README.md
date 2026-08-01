# Universal monotone selector search — July 29 continuation

This is the focused continuation workspace for the universal coalitionally
monotone core-selector problem on exact cooperative games.

**Status (2026-08-01): open.** Expanded searches over interlocking flows,
permutation-orbit couplings, the full signed exact-game tangent cone, and
multilevel lattices remained nonnegative. The best margins approach zero only
through a degenerate additive-game limit; no counterexample or universal proof
has been obtained.

## Orientation

- [`CURRENT_RESEARCH_REVIEW_2026-07-30.md`](CURRENT_RESEARCH_REVIEW_2026-07-30.md)
  is the best local status review.
- [`COMPLEMENTARY_FLOW_NO_GO.md`](COMPLEMENTARY_FLOW_NO_GO.md) records the main
  structural obstruction.
- [`CONNECTED_GRID7_SELF_SPLICE_RESULT.md`](CONNECTED_GRID7_SELF_SPLICE_RESULT.md)
  and [`N5_FULL_REACHABLE_BOUNDARY_RESULT.md`](N5_FULL_REACHABLE_BOUNDARY_RESULT.md)
  document the strongest exact families and boundary certificates.

The consolidated project, including broader selector counterexamples and the
main claim ledger, is
[`mfflanag1/exact-game-monotone-selection`](https://github.com/mfflanag1/exact-game-monotone-selection).

## Highest-value review

The main question is whether the protected-flow inequality suggested by the
finite exact certificates can be proved for arbitrary legal comparison
networks, or whether a genuinely interlocking topology escapes it. Additional
undirected brute-force search is lower value unless it changes the selector
basis or terminal incidence.

Treat numerical margins as search diagnostics. Any claimed counterexample
must be reconstructed over the rationals and checked for exact-game facets,
legal one-coordinate edges, and a negative compatibility certificate.
