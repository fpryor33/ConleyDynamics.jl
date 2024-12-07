
# Tip: Release Notes

Did you know you can add release notes too? Just add markdown formatted text underneath the comment after the text
"Release notes:" and it will be added to the registry PR, and if TagBot is installed it will also be added to the
release that TagBot creates. i.e.

@JuliaRegistrator register

Release notes:

## Breaking changes

- blah



## v0.1.1 (December 7, 2024)

- Removed unnecessary test from Lefschetz complex creation for speedup reasons.
- Minor changes to the documentation.

## v0.1.0 (November 29, 2024)

This is the initial official version hosted on the general Julia registry.

