# Connection Matrix Functions

## Connection Matrix Computation

```@docs
connection_matrix(::LefschetzComplex,::Vector{Vector{Int}};::Int,::Bool)
connection_matrix(::LefschetzComplex,::Vector{Vector{String}};::Int,::Bool)
cm_create!(::SparseMatrix,::Vector{Int};::Bool)
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

