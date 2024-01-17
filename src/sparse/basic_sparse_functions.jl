export sparse_size, sparse_low

"""
    n = sparse_size(matrix::SparseMatrix, dim::Int)

Number of rows (`dim=1`) or columns (`dim=2`) of a sparse matrix.
"""
function sparse_size(matrix::SparseMatrix, dim::Int)
    #
    # Return the number of columns or rows
    #

    if dim == 1
        return matrix.nrow
    elseif dim == 2
        return matrix.ncol
    else
        error("dim has to be 1 or 2!")
    end
end

"""
    n = sparse_low(matrix::SparseMatrix, col::Int)

Row index of the lowest nonzero matrix entry in column `col`.
"""
function sparse_low(matrix::SparseMatrix, col::Int)
    #
    # Return the location of the lowest nonzero entry in 
    # column col of the matrix
    #

    if length(matrix.columns[col]) == 0
        return 0
    else
        return last(matrix.columns[col])
    end
end

