export sparse_from_lists

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
function sparse_from_lists(nr::Int, nc::Int, typezero, typeone,
                           r::Vector{Int}, c::Vector{Int}, vals)
    #
    # Create a sparse matrix from a list of row indices, column indices,
    # as well as associated values. The values have to be in some field.
    #

    # Initialize the row, column, and entry vectors

    entries = Vector{Vector{typeof(vals[1])}}()
    columns = Vector{Vector{Int}}()
    rows = Vector{Vector{Int}}()

    for k=1:nc
        push!(entries,Vector{typeof(vals[1])}([]))
        push!(columns,Vector{Int}([]))
    end

    for k=1:nr
        push!(rows,Vector{Int}([]))
    end

    # Loop through the nonzero values and incorporate them
    
    nzindices = findall(x -> !(x==typezero), vals)

    for k in nzindices
        push!(entries[c[k]],vals[k])
        push!(columns[c[k]],r[k])
        push!(rows[r[k]],c[k])
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

    sm = SparseMatrix{typeof(vals[1])}(nr, nc, typezero, typeone,
                                       entries, columns, rows)

    # Return the struct

    return sm
end

