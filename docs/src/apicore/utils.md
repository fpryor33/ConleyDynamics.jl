# Utility Functions

```@meta
DocTestSetup = quote
    push!(LOAD_PATH,"../../../src/")
    using ConleyDynamics
end
```

## Simplicial Complexes

```@docs
create_simplicial_complex(::Vector{String},::Vector{Vector{Int}})
create_simplicial_complex(::Vector{String},::Vector{Vector{String}})
create_simplicial_rectangle(::Int,::Int)
create_cubical_complex(::Vector{String})
create_cubical_rectangle(::Int,::Int)
cube_field_size(::String)
cube_information(::String)
cube_label(::Int,::Int,::Vector{Int})
convert_simplices(::Vector{Vector{Int}},::Vector{String})
convert_simplices(::Vector{Vector{String}},::Vector{String})
```

## Lefschetz complexes

```@docs
lefschetz_closure(::LefschetzComplex,::Vector{Int})
lefschetz_closure(::LefschetzComplex,::Vector{String})
lefschetz_is_closed(::LefschetzComplex,::Vector{Int})
lefschetz_is_closed(::LefschetzComplex,::Vector{String})
lefschetz_clomo_pair(::LefschetzComplex,::Vector{Int})
lefschetz_clomo_pair(::LefschetzComplex,::Vector{String})
lefschetz_skeleton(::LefschetzComplex,::Vector{Int},::Int)
lefschetz_skeleton(::LefschetzComplex,::Vector{String},::Int)
lefschetz_skeleton(::LefschetzComplex,::Int)
lefschetz_subcomplex(::LefschetzComplex,::Vector{Int})
lefschetz_subcomplex(::LefschetzComplex,::Vector{String})
lefschetz_filtration(::LefschetzComplex,::Vector{Int})
lefschetz_filtration(::LefschetzComplex,::Vector{Vector{String}})
permute_lefschetz_complex(::LefschetzComplex,::Vector{Int})
```

## Multivector Fields

```@docs
convert_mvf(::Vector{Vector{Int}},::LefschetzComplex)
convert_mvf(::Vector{Vector{String}},::LefschetzComplex)
```

## Matrix Conversions

```@docs
convert_matrix_gfp(::SparseMatrix{Int},::Int)
```

