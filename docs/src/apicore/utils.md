# Utility Functions

```@meta
DocTestSetup = quote
    push!(LOAD_PATH,"../../../src/")
    using ConleyDynamics
end
```

## Simplicial Complexes

```@docs
create_simplicial_complex
create_simplicial_rectangle
create_simplicial_delaunay
```

## Cubical Complexes

```@docs
create_cubical_complex
create_cubical_rectangle
create_cubical_box
cube_field_size
cube_information
cube_label
```

## General Lefschetz complexes

```@docs
lefschetz_field
lefschetz_boundary
lefschetz_coboundary
lefschetz_openhull
lefschetz_closure
lefschetz_lchull
lefschetz_is_closed
lefschetz_is_locally_closed
lefschetz_clomo_pair
lefschetz_skeleton
lefschetz_subcomplex
lefschetz_closed_subcomplex
lefschetz_filtration
lefschetz_gfp_conversion
permute_lefschetz_complex
```

## Multivector Fields

```@docs
convert_cells
convert_cellsubsets
create_mvf_hull
create_planar_mvf
create_spatial_mvf
manifold_boundary
extract_multivectors
```

## General Helper Functions

```@docs
convert_planar_coordinates
convert_spatial_coordinates
planar_nontransverse_edges
```

