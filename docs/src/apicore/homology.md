# Homology Functions

## Regular Homology

```@docs
homology(::LefschetzComplex;::Int)
relative_homology(::LefschetzComplex,::Union{Vector{Int},Vector{String}};::Int)
```

## Persistent Homology

```@docs
persistent_homology(::LefschetzComplex,::Vector{Int};::Int)
```

## Reduction Algorithm

```@docs
ph_reduce!(::SparseMatrix;::Bool)
```

