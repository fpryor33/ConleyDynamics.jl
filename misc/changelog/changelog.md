
# Tip: Release Notes

Did you know you can add release notes too? Just add markdown formatted text
underneath the comment after the text "Release notes:" and it will be added to
the registry PR, and if TagBot is installed it will also be added to the release
that TagBot creates. i.e., you could add a "## Breaking Changes" header with a
description, etc.

*******************************************************************************

@JuliaRegistrator register

Release notes:

## v0.1.6 (March 1, 2025)

- Added triangulations for the torus, the Klein bottle and the
  projective plane.
- Updated the documentation.

## v0.1.5 (February 25, 2025)

- Added `create_random_filter` to create a random injective filter.
- Added `filter_shallow_pairs` to find all shallow pairs of a filter.
- Added `filter_induced_mvf` to compute the multivector field induced
  by a filter.

## v0.1.4 (February 21, 2025)

- Added `lefschetz_reduction` to construct a smaller Lefschetz complex
  via elementary reductions.
- Added `sparse_get_nz_row` to extract all column indices corresponding
  to nonzero entries in a row.
- Added `sparse_is_zero` to check whether a sparse matrix is zero.
- Added consistency check to `create_lefschetz_gf2`.

## v0.1.3 (February 17, 2025)

- Added `create_lefschetz_gf2` function for quick Lefschetz complex creation
  over a field of characteristic 2.
- Simplified some examples in the documentation based on the new function.

## v0.1.2 (February 12, 2025)

- Fixed a rare exception in the multivector field creation

## v0.1.1 (December 7, 2024)

- Removed unnecessary test from Lefschetz complex creation for speedup reasons.
- Minor changes to the documentation.

## v0.1.0 (November 29, 2024)

This is the initial official version hosted on the general Julia registry.

