export sparse_from_full, full_from_sparse

"""
    sm = sparse_from_full(matrix::Matrix{Int})

Create sparse matrix from full matrix.
"""
function sparse_from_full(matrix::Matrix{Int})
    #
    # Create sparse matrix from full matrix

    # Initialize the row, column, and entry vectors

    nr = size(matrix,1)
    nc = size(matrix,2)
    tzero = Int(0)
    tone  = Int(1)
    
    # Create the lists for the nonzero entries

    r = Vector{Int}([])
    c = Vector{Int}([])
    vals = Vector{Int}([])

    for k=1:nr
        for m=1:nc
            mvalue = matrix[k,m]
            if !(mvalue == 0)
                push!(r,k)
                push!(c,m)
                push!(vals,mvalue)
            end
        end
    end

    # Create the sparse matrix and return it

    sm = sparse_from_lists(nr, nc, tzero, tone, r, c, vals)
    return sm
end

"""
    fm = full_from_sparse(sm::SparseMatrix)

Create full matrix from sparse matrix.
"""
function full_from_sparse(sm::SparseMatrix)
    #
    # Create full matrix from sparse matrix
    #

    # Initialize the full matrix with zeros

    fm = fill(sm.zero, (sm.nrow, sm.ncol))

    # Extract the nonzero elements

    nr, nc, tzero, tone, r, c, vals = lists_from_sparse(sm)

    # Set the nonzero elements

    for k=1:length(r)
        fm[r[k],c[k]] = vals[k]
    end

    # Return the full matrix
    
    return fm
end

