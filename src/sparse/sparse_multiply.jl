export sparse_multiply

"""
    sparse_multiply(A::SparseMatrix, B::SparseMatrix)

Multiply two sparse matrices.

Exceptions are raised if the matrix product is not defined 
or the entry types do not match.
"""
function sparse_multiply(A::SparseMatrix, B::SparseMatrix)
    #
    # Multiply two sparse matrices.
    #

    # Check whether the product exists

    if !(A.ncol == B.nrow)
        error("Matrix product undefined!")
    end

    if !(typeof(A.zero) == typeof(B.zero))
        error("The sparse matrices need to have the same type!")
    end

    if !(typeof(A.char) == typeof(B.char))
        error("The sparse matrices have to be over the same field!")
    end

    # Perform the matrix product

    r = Vector{Int}()
    c = Vector{Int}()
    v = Vector{typeof(A.zero)}()
    tzero = A.zero

    for k=1:A.nrow
        for m=1:B.ncol
            dpindex = intersect(A.rows[k],B.columns[m])
            if length(dpindex) > 0
                dotprod = tzero
                for j in dpindex
                    dotprod += A[k,j] * B[j,m]
                end
                if !(dotprod == tzero)
                    push!(r,k)
                    push!(c,m)
                    push!(v,dotprod)
                end
            end
        end
    end

    # Construct and return the sparse product matrix

    sm = sparse_from_lists(A.nrow, B.ncol, A.char, A.zero, A.one, r, c, v)
    return sm
end

"""
    Base.:*(A::SparseMatrix, B::SparseMatrix)

Multiply two sparse matrices.

Exceptions are raised if the matrix product is not defined 
or the entry types do not match.
"""
Base.:*(A::SparseMatrix,B::SparseMatrix) = sparse_multiply(A,B)

