export update_low!

"""
    update_low!(matrix, lowvec; startindex=1)

Compute the indices of the lowest nonzero entries in each column of `matrix`.

The indices are stored in `lowvec`. If `startindex > 1`, then only the values
starting with column `startindex` are updated.
"""
function update_low!(matrix, lowvec; startindex=1)
    #
    # Compute the indices of the lowest nonzero
    # entries in each column of the matrix. The
    # indices are stored in lowvec. If one chooses
    # startindex > 1, then only the values starting
    # with column startindex are updated.
    #
    numcolumns    = size(matrix, 2)
    numrows       = size(matrix, 1)
    
    for j = startindex:numcolumns
        foundit = false
        k = numrows
        while (!foundit) && (k >= 1)
            if matrix[k,j] != 0
                foundit = true
            else
                k -= 1
            end
        end
        lowvec[j] = k
    end
    
    return nothing
end

