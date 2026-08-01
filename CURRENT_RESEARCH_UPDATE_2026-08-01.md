# Research update: universal coalitionally monotone core selection

**Updated:** 2026-08-01  
**Central status:** open

## Executive result

The central universal-selector question is not yet solved, but the search
frontier moved decisively. Two exact results now eliminate every protected
block-flow pseudoforest and the entire most dangerous common-grand
five-player \(K_{2,3}\) class:

1. Every weighted alternating protected cycle is safe in arbitrary
   dimension, even with unequal grand worths. Together with leaf contraction,
   this closes every forest or unicyclic component.
2. All 20 irreducible maximally sharp five-player \(K_{2,3}\) symmetry orbits
   have a common core point whenever the grand worth is common. The
   classification covers 4,440 sink-normalized configurations and 66,108
   exact rational envelope certificates, with zero failures.

Thus the first cycle-rank-two class survives neither local premium alignment
nor global protected coupling. Remaining candidates must use unequal grand
worths, non-sharp or differently oriented rank-two components, higher cycle
rank, or more players.

## Basic explanation

A core allocation pays every coalition at least what that coalition can earn
alone. An exact game is especially rich: every coalition is paid exactly its
worth at some core allocation.

The desired rule must choose one core allocation from every exact game. If a
single coalition's worth rises, none of its members may be paid less. A finite
failure would appear as a network of games whose protected player flows have
positive total lower expectation.

Trees and single cycles cannot do this: exactness lets coalition-worth terms
be paired along the flow, where they cancel. The next graph is \(K_{2,3}\),
which has two independent cycles. The locally sharp five-player versions were
designed so every terminal has the largest possible non-additivity premium.
Nevertheless, the protected equalities force the three upper games to have a
balanced pointwise envelope. Their cores therefore intersect, and that common
allocation makes every zero-sum flow objective nonpositive.

## Exact outputs

- [Alternating-cycle theorem](ALTERNATING_CYCLE_NO_GO_RESULT.md)
- [Sharp K2,3 common-core classification](N5_K23_COMMON_CORE_RESULT.md)
- [Mixed-divergence premium classification](N5_MIXED_DIVERGENCE_PREMIUM_RESULT.md)
- [Johnson-square exact no-go](N5_JOHNSON_SQUARE_NO_GO_RESULT.md)
- [K2,3 classification summary](results/n5_k23_all_sharp_common_core_classification_exact.json)
- [Exact balanced-vertex catalogue](results/n5_balanced_vertices_exact.json)
- [Exact empty-common-core negative control](results/n5_k23_empty_common_core_witness_exact.json)

The sharp \(K_{2,3}\) verification chain is:

- 4,440 configurations after fixing the first row's sink labels;
- 20 orbits after player, source, and sink symmetries;
- 1,291 extreme five-player balanced collections;
- 66,108 envelope assignments across all orbits;
- 66,108 rational LP dual reconstructions;
- zero reconstruction failures; and
- exact common-core gap zero in all 20 orbits.

The exact local premium archives were rechecked independently in the same
run: 2,435 three-level branches, 2,296 four-level branches, and 513 branches
of the zero-premium 2-1-2 identity all pass.

## Negative control and methodological limit

The common-core method is not universal. A non-sharp \(K_{2,3}\) with masks

\[
\begin{pmatrix}4&15&8\\4&23&8\end{pmatrix}
\]

admits an exact unanimity-game realization whose common-core budget gap is
one. Yet its terminal selector objective is \(-3\), and direct optimization
returns zero as the topology maximum. This cleanly separates two notions:
having no constant allocation is much weaker than having no monotone
selection. The common-core separation also has the wrong flow orientation:
its positive terms occur at upper sink games.

## Verification depth and level of work

The new results use exact integer enumeration, the complete 280-facet
five-player exact-game cone at every node, exact protected-path relations, and
rational dual stationarity/sign/objective checks. They have extensive
internal and independent-code verification but no external peer review.

- Thesis chapter: exceeded.
- PhD thesis scale: yes, for the combined project.
- Professional-research level: yes; the dimension-free cycle theorem and the
  complete \(K_{2,3}\) classification are research-grade results.
- Breakthrough: not yet on the central universal question. The new theorem is
  a substantial structural advance, but it closes a sharply defined class
  rather than the universal quantifier.

## Strongest next steps

1. Prove or disprove that unequal grand worths can escape the sharp
   \(K_{2,3}\) common-grand theorem. Numerical sweeps of all 20 orbits remain
   zero, but this is not exact.
2. Characterize which non-sharp \(K_{2,3}\) empty-core separations can align
   with forward protected flow. The first exact negative control cannot.
3. Extend the balanced-envelope method to other rank-two orientations, then
   to \(K_{3,3}\) or six-player degree-three terminals.
4. If any positive floating obstruction appears, immediately reconstruct the
   games and flow over the rationals and lift every protected relation to a
   legal one-coordinate exact-game path before claiming a counterexample.
5. In parallel, search for a symbolic charging theorem that converts every
   higher-cycle local premium into protected-path slack. Such a theorem would
   be the clearest route toward universal existence.
