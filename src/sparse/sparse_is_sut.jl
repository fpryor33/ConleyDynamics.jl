export sparse_is_sut

"""
    bool = sparse_is_sut(sm::SparseMatrix)

Check whether the sparse matrix is strictly upper triangular.
"""
function sparse_is_sut(sm::SparseMatrix)
    #
    # Is the matrix strictly upper triangular?
    #

    for k=1:sm.nrow
        if length(sm.rows[k]) > 0
            if sm.rows[k][1] <= k
                return false
            end
        end
    end

    return true
end

