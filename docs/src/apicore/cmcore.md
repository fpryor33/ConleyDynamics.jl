# Connection Matrix Functions

## Connection Matrix Computation

```@docs
connection_matrix(::LefschetzComplex,::Vector{Vector{Int}};::Int,::Bool)
connection_matrix(::LefschetzComplex,::Vector{Vector{String}};::Int,::Bool)
cm_create!(::Matrix,::Vector{Int};::Bool)
cm_create!(::SparseMatrix,::Vector{Int};::Bool)
```

## Poset Order Functions

```@docs
admissible_order(bndmatrix::Matrix{Int}, mvf::Vector{Vector{Int}})
admissible_order(bndmatrix::SparseMatrix{Int}, mvf::Vector{Vector{Int}})
renumber_poset!(poset::Vector{Int})
```

## Sparse Matrix Helper Functions

```@docs
cm_columns(matrix::SparseMatrix, psetvec::Vector{Int})
homogeneous_columns(::SparseMatrix,::Vector{Int})
is_homogeneous(::SparseMatrix,::Vector{Int},::Int)
target_columns(matrix::SparseMatrix, psetvec::Vector{Int})
```

## Full Matrix Helper Functions

```@docs
cm_columns(matrix::Matrix, lowvec::Vector{Int}, psetvec::Vector{Int})
homogeneous_columns(::Matrix,::Vector{Int},::Vector{Int})
is_homogeneous(::Matrix,::Vector{Int},::Vector{Int},::Int)
target_columns(::Matrix,::Vector{Int},::Vector{Int})
update_low!(matrix::Matrix, lowvec::Vector{Int}; startindex::Int=1)
```

