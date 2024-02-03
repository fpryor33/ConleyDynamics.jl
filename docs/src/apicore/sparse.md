# Sparse Matrix Functions

## Access Functions

```@docs
sparse_get_entry(matrix::SparseMatrix, ri::Int, ci::Int)
Base.getindex(matrix::SparseMatrix, ri::Int, ci::Int)
sparse_set_entry!(matrix::SparseMatrix, ri::Int, ci::Int, val)
Base.setindex!(matrix::SparseMatrix, val, ri::Int, ci::Int)
sparse_get_column(matrix::SparseMatrix, ci::Int)
sparse_get_nz_column(matrix::SparseMatrix, ci::Int)
sparse_minor(sm::SparseMatrix, rvec::Vector{Int}, cvec::Vector{Int})
```

## Basic Functions

```@docs
sparse_size(matrix::SparseMatrix, dim::Int)
sparse_low(matrix::SparseMatrix, col::Int)
sparse_is_sut(sm::SparseMatrix)
sparse_identity(n::Int, tone)
sparse_fullness(sm::SparseMatrix)
sparse_sparsity(sm::SparseMatrix)
sparse_show(sm::SparseMatrix)
sparse_show(sm::SparseMatrix{Int})
Base.show(io::IO, sm::SparseMatrix)
Base.show(io::IO, sm::SparseMatrix{Int})
```

## Elementary Matrix Operations

```@docs
sparse_add_column!(matrix::SparseMatrix, ci1::Int, ci2::Int, c)
sparse_add_row!(matrix::SparseMatrix, ri1::Int, ri2::Int, c)
sparse_permute(sm::SparseMatrix, pr::Vector{Int}, pc::Vector{Int})
sparse_remove!(matrix::SparseMatrix, ri::Int, ci::Int)
sparse_multiply(A::SparseMatrix,B::SparseMatrix)
Base.:*(::SparseMatrix,::SparseMatrix)
```

## Conversion Functions

```@docs
sparse_from_full(matrix::Matrix{Int})
full_from_sparse(sm::SparseMatrix)
sparse_from_lists(nr::Int, nc::Int, tzero, tone, r::Vector{Int}, c::Vector{Int}, vals)
lists_from_sparse(sm::SparseMatrix)
```

