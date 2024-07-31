# Lefschetz Complex Functions

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
get_cubical_coords
```

## Lefschetz Complex Queries

```@docs
lefschetz_field
lefschetz_is_closed
lefschetz_is_locally_closed
```

## Topological Features

```@docs
lefschetz_boundary
lefschetz_coboundary
lefschetz_closure
lefschetz_openhull
lefschetz_lchull
lefschetz_clomo_pair
lefschetz_skeleton
manifold_boundary
```

## Lefschetz Subcomplexes

```@docs
lefschetz_subcomplex
lefschetz_closed_subcomplex
permute_lefschetz_complex
```

## Lefschetz Helper Functions

```@docs
lefschetz_gfp_conversion
lefschetz_filtration
```

## Cell Subset Helper Functions

```@docs
convert_cells
convert_cellsubsets
```

## Coordinate Helper Functions

```@docs
convert_planar_coordinates
convert_spatial_coordinates
```

