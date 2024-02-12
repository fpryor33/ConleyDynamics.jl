export sparse_from_lists, lists_from_sparse

"""
    sparse_from_lists(nr, nc, tchar, tzero, tone, r, c, v)

Create sparse matrix from lists describing the entries.

The vectors `r`, `c`, and `v` have to have the same length
and the matric has entry `v[k]` at `(r[k],c[k])`. Zero
entries will be ignored, and multiple entries for the same
matrix position raise an error.

The input arguments have the following meaning:
* `nr::Int`: Number of rows
* `nc::Int`: Number of columns
* `tchar`: Field characteristic if `T==Int`
* `tzero::T`: Number 0 of type `T`
* `tone::T`:  Number 1 of type `T`
* `r::Vector{Int}`: Vector of row indices
* `c::Vector{Int}`: Vector of column indices
* `v::Vector{T}`: Vector of matrix entries
If `tchar>0`, then the entries in `v` are all
replaced by their values mod `tchar`.
"""
function sparse_from_lists(nr::Int, nc::Int, tchar::Int, tzero, tone,
                           r::Vector{Int}, c::Vector{Int}, v)
    #
    # Create a sparse matrix from a list of row indices, column indices,
    # as well as associated values. The values have to be in some field.
    #

    # Make sure that we are in GF(p) if tchar > 0

    if (tzero isa Int) && (tchar > 0)
        vals = Vector{Int}()
        for k=1:length(v)
            push!(vals,mod(v[k], tchar))
        end
    else
        vals = v
    end

    # Initialize the row, column, and entry vectors

    entries = [Vector{typeof(tzero)}([]) for _ in 1:nc]
    columns = [Vector{Int}([]) for _ in 1:nc]
    rows    = [Vector{Int}([]) for _ in 1:nr]

    # If the lists have length zero, return the zero matrix

    if length(r)==0
        sm = SparseMatrix{typeof(tzero)}(nr, nc, tchar, tzero, tone,
                                         entries, columns, rows)
        return sm
    end

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

    # Create and return the sparse matrix object

    sm = SparseMatrix{typeof(tzero)}(nr, nc, tchar, tzero, tone,
                                     entries, columns, rows)
    return sm
end

"""
    nr, nc, tchar, tzero, tone, r, c, v = lists_from_sparse(sm::SparseMatrix)

Create list representation from sparse matrix.

The output variables are exactly what is needed to create
a sparse matrix object using `sparse_from_lists`.
"""
function lists_from_sparse(sm::SparseMatrix)
    #
    # Create list representation from sparse matrix
    #

    # Extract the constant fields of the struct

    nr = sm.nrow
    nc = sm.ncol
    tchar = sm.char
    tzero = sm.zero
    tone  = sm.one

    # Create the lists

    r = Vector{Int}([])
    c = Vector{Int}([])
    v = Vector{typeof(tzero)}([])

    for k=1:nc
        for m=1:length(sm.columns[k])
            push!(r,sm.columns[k][m])
            push!(c,k)
            push!(v,sm.entries[k][m])
        end
    end

    # Return the results

    return nr, nc, tchar, tzero, tone, r, c, v
end

