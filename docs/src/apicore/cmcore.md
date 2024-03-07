# Conley Theory Functions

## Conley Index Computations

```@docs
conley_index(::LefschetzComplex,::Vector{Int};::Int)
conley_index(::LefschetzComplex,::Vector{String};::Int)
```

## Connection Matrix Computation

```@docs
connection_matrix(::LefschetzComplex,::MultiVectorField;::Int,::Bool)
cm_reduce!(::SparseMatrix,::Vector{Int};::Bool)
```

## Poset Order Functions

```@docs
admissible_order(bndmatrix::SparseMatrix, mvf::Vector{Vector{Int}})
renumber_poset!(poset::Vector{Int})
```

## Matrix Helper Functions

```@docs
cm_columns(matrix::SparseMatrix, psetvec::Vector{Int})
homogeneous_columns(::SparseMatrix,::Vector{Int})
is_homogeneous(::SparseMatrix,::Vector{Int},::Int)
target_columns(matrix::SparseMatrix, psetvec::Vector{Int})
```

