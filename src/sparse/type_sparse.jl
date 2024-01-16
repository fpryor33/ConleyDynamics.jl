export SparseMatrix

"""
    SparseMatrix{T}

Composite data type for a sparse matrix with entries of type `T`.

The struct has the following fields:
* `const nrow::Int`:
* `const ncol::Int`:
* `const zero::T`:
* `const one::T`:
* `entries::Vector{Vector{T}}`:
* `columns::Vector{Vector{Int}}`:
* `rows::Vector{Vector{Int}}`:
"""
mutable struct SparseMatrix{T}
    const nrow::Int
    const ncol::Int
    const zero::T
    const one::T
    entries::Vector{Vector{T}}
    columns::Vector{Vector{Int}}
    rows::Vector{Vector{Int}}
end

