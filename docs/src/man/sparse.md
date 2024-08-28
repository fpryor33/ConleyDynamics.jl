# Sparse Matrices

While Julia provides a data structure for sparse matrix computations, the
employed design decisions make it difficult to use this implementation for
computations over finite fields. This is mainly due to the fact that in
the Julia implementation, it is assumed that one can determine the zero and
one elements from the data type alone. However, a finite field data type
generally also depends on additional parameters, such as the characteristic
of the field.

Since the algorithms underlying `ConleyDynamics.jl` only require basic row
and column operations, a specialized sparse matrix implementation is provided
in the package. It is briefly described in the following.

## Sparse Matrix Format

```@docs; canonical=false
SparseMatrix
```

## Creating Sparse Matrices

- [`sparse_from_full`](@ref)
- [`full_from_sparse`](@ref)
- [`sparse_from_lists`](@ref)
- [`lists_from_sparse`](@ref)
- [`sparse_identity`](@ref)

## Sparse Matrix Access

- [`sparse_get_entry`](@ref)
- [`sparse_set_entry!`](@ref)
- [`sparse_get_column`](@ref)
- [`sparse_get_nz_column`](@ref)
- [`sparse_minor`](@ref)

You can also use `y = A[i,j]` and `A[i,j] = val`.

## Sparse Matrix Information

- [`sparse_size`](@ref)
- [`sparse_low`](@ref)
- [`sparse_is_sut`](@ref)
- [`sparse_fullness`](@ref)
- [`sparse_sparsity`](@ref)
- [`sparse_show`](@ref)

## Elementary Matrix Operations

- [`sparse_add_column!`](@ref)
- [`sparse_add_row!`](@ref)
- [`sparse_permute`](@ref)
- [`sparse_remove!`](@ref)
- [`sparse_multiply`](@ref)

You can also use `A*A` to compute the product of sparse matrices.






