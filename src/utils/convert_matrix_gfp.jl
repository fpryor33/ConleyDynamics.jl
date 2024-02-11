export convert_matrix_gfp

"""
    convert_matrix_gfp(matrix::SparseMatrix{Int}, p::Int)

Convert a sparse integer matrix to a finite field sparse matrix
over `GF(p)`. For the choice `p=0` the rationals are used.
"""
function convert_matrix_gfp(matrix::SparseMatrix{Int}, p::Int)
    #
    # Convert an integer matrix to a finite field matrix over GF(p).
    #

    # If the sparse matrix already has p>0, either raise an
    # error or pass the matrix through

    if matrix.char > 0
        if matrix.char == p
            return matrix
        else
            error("Cannot convert matrix from GF(q) to GF(p)!")
        end
    end

    # Convert sparse matrix to lists

    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(matrix)

    # Convert zero, one, and the matrix entries to GF(p)

    if (p > 0)
        gfpchar = p
        gfpzero = Int(0)
        gfpone  = Int(1)
        gfpvals = Vector{Int}()
        for k=1:length(vals)
            push!(gfpvals,vals[k]) # mod will be done by sparse_from_lists
        end
    elseif (p == 0)
        gfpchar = 0
        gfpzero = Rational{Int}(0)
        gfpone  = Rational{Int}(1)
        gfpvals = Vector{Rational{Int}}()
        for k=1:length(vals)
            push!(gfpvals,Rational{Int}(vals[k]))
        end
    else
        error("Wrong characteristic p!")
    end

    # Create and return the sparse GFP(p) matrix

    gfpmatrix = sparse_from_lists(nr, nc, gfpchar, gfpzero, gfpone, r, c, gfpvals)
    return gfpmatrix
end

