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

Sparse matrices in this package have to be of the composite
data type [`SparseMatrix`](@ref), which is structured as follows:

```@docs; canonical=false
SparseMatrix
```

In this struct, the type `T` has to be either `Int` or `Rational{Int}`,
depending on whether the sparse matrix is interpreted as a matrix
with entries in the finite field ``GF(p)`` for some prime ``p``, or
over the field of rationals, respectively. The data type has the
following fields:

- `nrow::Int` designates the number of rows.
- `ncol::Int` gives the number of columns.
- `char::Int` specifies the characteristic of
  the underlying field ``F``. If `char=0`, then the field is the
  rationals ``\mathbb{Q}``, and one has to have `T = Rational{Int}`.
  If, on the other hand, the finite field ``F = GF(p)`` is used,
  then `char=p` has to be a prime number. In this case, the data
  type of the matrix entries has to be `T = Int`.
- `zero::T` provides 0 in the data type `T`.
- `one::T` provides 1 in the data type `T`.
- `columns::Vector{Vector{Int}}` is a vector of integer vectors, which
  contains the row indices of nonzero matrix entries in each column.
  More precisely, `columns[k]` contains an increasing list of row indices,
  which give the locations of all nonzero entries in column `k`. Note that
  the list for each colum has to be strictly increasing.
- `rows::Vector{Vector{Int}}` is a vector of integer vectors, which
  contains the column indices of nonzero matrix entries in each row.
  It is the precise dual to the previous field. This time, `rows[k]`
  contains an increasing list of column indices, which correspond to the
  nonzero entries of the matrix in the `k`-th row.
- `entries::Vector{Vector{T}}` is a vector of vectors which contains 
  the actual matrix entries. It is organized in exactly the same way as
  the field `columns`. In other words, for every `k = 1,..,ncol` the
  matrix entry in column `k` and row `columns[k][j]` is given by
  `entries[k][j]`, where `j` indexes the nonzero column entries from
  top to bottom.

This data structure is clearly redundant, in the sense that the
field `rows` is not needed to uniquely determine the matrix. However,
the type [`SparseMatrix`](@ref) is fundamental for almost every aspect
of `ConleyDynamics.jl`, as it is used to encode the incidence coefficient
map ``\kappa``, and therefore also the matrix representation of the
boundary operator ``\partial``. And for many operations on or queries
of Lefschetz complexes, one needs fast access to both the cells in
the boundary and the coboundary of a given cells. While the boundary
can easily be accessed via the field `columns`, the fast coboundary
access is aided by the field `rows`.

We would like to point out that in view of the different underlying
fields, sparse matrices should only be manipulated using the specific
commands provided by the package. These are described in detail below.
If there is a need for additional functionality beyond these first
methods, it can be added at a later point in time.

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

## Elementary Matrix Operations

- [`sparse_add_column!`](@ref)
- [`sparse_add_row!`](@ref)
- [`sparse_permute`](@ref)
- [`sparse_remove!`](@ref)
- [`sparse_multiply`](@ref)

You can also use `A*A` to compute the product of sparse matrices.

## Sparse Matrix Information

- [`sparse_size`](@ref)
- [`sparse_low`](@ref)
- [`sparse_is_sut`](@ref)
- [`sparse_fullness`](@ref)
- [`sparse_sparsity`](@ref)
- [`sparse_show`](@ref)

