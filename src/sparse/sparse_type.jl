export SparseMatrix

"""
    SparseMatrix{T}

Composite data type for a sparse matrix with entries of type `T`.

The struct has the following fields:
* `const nrow::Int`: Number of rows
* `const ncol::Int`: Number of columns
* `const char::Int`: Characteristic of type `T`
* `const zero::T`: Number 0 of type `T`
* `const one::T`:  Number 1 of type `T`
* `entries::Vector{Vector{T}}`: Matrix entries corresponding to `columns`
* `columns::Vector{Vector{Int}}`: `columns[k]` points to nonzero entries in column k
* `rows::Vector{Vector{Int}}`: `rows[k]` points to nonzero entries in the k-th row
"""
mutable struct SparseMatrix{T}
    const nrow::Int
    const ncol::Int
    const char::Int
    const zero::T
    const one::T
    entries::Vector{Vector{T}}
    columns::Vector{Vector{Int}}
    rows::Vector{Vector{Int}}
end


# Specialized deepcopy for sparse matrices. The generic deepcopy records
# every inner vector in an IdDict, which is slow for large matrices. The
# inner vectors of a sparse matrix are never shared, so it suffices to
# copy them one by one. Entries of non-isbits type are deep copied.

function Base.deepcopy_internal(sm::SparseMatrix{T}, stackdict::IdDict) where T
    haskey(stackdict, sm) && return stackdict[sm]::SparseMatrix{T}
    if isbitstype(T)
        entries = [copy(e) for e in sm.entries]
    else
        entries = Base.deepcopy_internal(sm.entries, stackdict)
    end
    columns = [copy(c) for c in sm.columns]
    rows    = [copy(r) for r in sm.rows]
    smcopy  = SparseMatrix{T}(sm.nrow, sm.ncol, sm.char, sm.zero, sm.one,
                              entries, columns, rows)
    stackdict[sm] = smcopy
    return smcopy
end
