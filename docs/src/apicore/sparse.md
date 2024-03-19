# Sparse Matrix Functions

## Internal Sparse Matrix Representation

```@docs
SparseMatrix
```

## Access Functions

```@docs
sparse_get_entry
Base.getindex(matrix::SparseMatrix, ri::Int, ci::Int)
sparse_set_entry!
Base.setindex!(matrix::SparseMatrix, val, ri::Int, ci::Int)
sparse_get_column
sparse_get_nz_column
sparse_minor
```

## Basic Functions

```@docs
sparse_size
sparse_low
sparse_is_sut
sparse_identity
sparse_fullness
sparse_sparsity
sparse_show
```

## Elementary Matrix Operations

```@docs
sparse_add_column!
sparse_add_row!
sparse_permute
sparse_remove!
sparse_multiply
Base.:*(::SparseMatrix,::SparseMatrix)
```

## Conversion Functions

```@docs
sparse_from_full
full_from_sparse
sparse_from_lists
lists_from_sparse
```

