# Utility Functions

## Simplicial Complexes

```@docs
create_simplicial_complex(::Vector{String},::Vector{Vector{Int}})
create_simplicial_complex(::Vector{String},::Vector{Vector{String}})
convert_simplices(::Vector{Vector{Int}},::Vector{String})
convert_simplices(::Vector{Vector{String}},::Vector{String})
```

## Lefschetz complexes

```@docs
lefschetz_closure(::LefschetzComplex,::Vector{Int})
lefschetz_closure(::LefschetzComplex,::Vector{String})
permute_lefschetz_complex(::LefschetzComplex,::Vector{Int})
convert_lefschetz_sparse(::LefschetzComplex)
```

## Multivector Fields

```@docs
convert_mvf(::Vector{Vector{Int}},::LefschetzComplex)
convert_mvf(::Vector{Vector{String}},::LefschetzComplex)
```

## Matrix Conversions

```@docs
convert_matrix_gfp(::Matrix{Int},::Int)
convert_matrix_gfp(::SparseMatrix{Int},::Int)
convert_matrix_int(::Matrix)
```

