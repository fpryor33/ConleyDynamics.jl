export homogeneous_columns, is_homogeneous

"""
    homogeneous_columns(matrix::SparseMatrix, psetvec::Vector{Int})
 
Determine which columns of `matrix` are homogenenous columns.
"""
function homogeneous_columns(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Determine which columns are homogenous columns
    #
    numcolumns = sparse_size(matrix, 2)
    hcols      = Vector(zeros(Bool,numcolumns))

    for j = 1:numcolumns
        hcols[j] = is_homogeneous(matrix, psetvec, j)
    end

    return hcols
end

"""
    is_homogeneous(matrix::SparseMatrix, psetvec::Vector{Int}, cindex::Int)

Decide whether a column of `matrix` is homogeneous.
"""
function is_homogeneous(matrix::SparseMatrix, psetvec::Vector{Int},
                        cindex::Int)
    #
    # Decide whether a column is homogeneous.
    #

    lowc = sparse_low(matrix,cindex)

    if lowc == 0
        return false
    else
        targetindex = lowc
        if psetvec[cindex] == psetvec[targetindex]
            return true
        else
            return false
        end
    end
end

