export sparse_from_lists

"""
    sm = sparse_from_lists(nr, nc, tzero, tone, r, c, vals)

Create sparse matrix from lists describing the entries.

The vectors `r`, `c`, and `vals` have to have the same length
and the matric has entry `vals[k]` at `(r[k],c[k])`. Zero
entries will be ignored, and multiple entries for the same
matrix position raise an error.

The input arguments have the following meaning:
* `nr::Int`: Number of rows
* `nc::Int`: Number of columns
* `tzero::T`: Number 0 of type `T`
* `tone::T`:  Number 1 of type `T`
* `r::Vector{Int}`: Vector of row indices
* `c::Vector{Int}`: Vector of column indices
* `vals::Vector{T}`: Vector of matrix entries
"""
function sparse_from_lists(nr::Int, nc::Int, tzero, tone,
                           r::Vector{Int}, c::Vector{Int}, vals)
    #
    # Create a sparse matrix from a list of row indices, column indices,
    # as well as associated values. The values have to be in some field.
    #

    # Initialize the row, column, and entry vectors

    entries = [Vector{typeof(vals[1])}([]) for _ in 1:nc]
    columns = [Vector{Int}([]) for _ in 1:nc]
    rows    = [Vector{Int}([]) for _ in 1:nr]

    # Loop through the nonzero values and incorporate them
    
    nzindices = findall(x -> !(x==tzero), vals)

    for k in nzindices
        push!(entries[c[k]],vals[k])
        push!(columns[c[k]],r[k])
        push!(rows[r[k]],c[k])
    end

    # Check whether the entries are all unique

    for k=1:length(columns)
        if length(columns[k]) > length(unique(columns[k]))
            error("Duplicate sparse matrix entries!")
        end
    end

    # Sort the representations

    for k=1:length(columns)
        sp = sortperm(columns[k])
        columns[k] = columns[k][sp]
        entries[k] = entries[k][sp]
    end

    for k=1:length(rows)
        sort!(rows[k])
    end

    # Create the sparse matrix object

    sm = SparseMatrix{typeof(vals[1])}(nr, nc, tzero, tone,
                                       entries, columns, rows)

    # Return the struct

    return sm
end

