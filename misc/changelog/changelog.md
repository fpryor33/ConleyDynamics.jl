
# Tip: Release Notes

Did you know you can add release notes too? Just add markdown formatted text
underneath the comment after the text "Release notes:" and it will be added to
the registry PR, and if TagBot is installed it will also be added to the release
that TagBot creates. i.e., you could add a "## Breaking Changes" header with a
description, etc.

*******************************************************************************

@JuliaRegistrator register

Release notes:

## v0.1.4 ()

- Added `sparse_is_zero` to check whether a sparse matrix is zero.

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

