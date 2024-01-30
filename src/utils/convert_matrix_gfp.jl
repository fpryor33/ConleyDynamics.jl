export convert_matrix_gfp

"""
    gfpmatrix = convert_matrix_gfp(matrix::Matrix{Int}, p::Int)

Convert an integer matrix to a finite field matrix over GF(p).
"""
function convert_matrix_gfp(matrix::Matrix{Int}, p::Int)
    #
    # Convert an integer matrix to a finite field matrix over GF(p).
    #
    m,n = size(matrix)
    FF = GF(p)
    gfpmatrix = [FF(matrix[i,j]) for i=1:m, j=1:n]
    return gfpmatrix
end

"""
    gfpmatrix = convert_matrix_gfp(matrix::SparseMatrix{Int}, p::Int)

Convert a sparse integer matrix to a finite field sparse matrix over GF(p).
"""
function convert_matrix_gfp(matrix::SparseMatrix{Int}, p::Int)
    #
    # Convert an integer matrix to a finite field matrix over GF(p).
    #

    # Convert sparse matrix to lists

    nr, nc, tzero, tone, r, c, vals = lists_from_sparse(matrix)

    # Convert zero, one, and the matrix entries to GF(p)

    FF = GF(p)

    gfpzero = FF(tzero)
    gfpone  = FF(tone)

    gfpvals = Vector{typeof(gfpzero)}()
    for k=1:length(vals)
        push!(gfpvals,FF(vals[k]))
    end

    # Create and return the sparse GFP(p) matrix

    gfpmatrix = sparse_from_lists(nr, nc, gfpzero, gfpone, r, c, gfpvals)
    return gfpmatrix
end

