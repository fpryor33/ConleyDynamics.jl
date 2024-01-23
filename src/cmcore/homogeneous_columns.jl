export homogeneous_columns, is_homogeneous

"""
    hcols = homogeneous_columns(matrix, lowvec, psetvec)
 
Determine which columns of `matrix` are homogenenous columns.
"""
function homogeneous_columns(matrix::Matrix, lowvec::Vector{Int},
                             psetvec::Vector{Int})
    #
    # Determine which columns are homogenous columns
    #
    numcolumns = size(matrix, 2)
    hcols      = Vector(zeros(Bool,numcolumns))

    for j = 1:numcolumns
        hcols[j] = is_homogeneous(matrix, lowvec, psetvec, j)
    end

    return hcols
end

"""
    hcols = homogeneous_columns(matrix::SparseMatrix, psetvec)
 
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
    bool = is_homogeneous(matrix, lowvec, psetvec, cindex)

Decide whether a column of `matrix` is homogeneous.
"""
function is_homogeneous(matrix::Matrix, lowvec::Vector{Int},
                        psetvec::Vector{Int}, cindex::Int)
    #
    # Decide whether a column is homogeneous. The function expects that the
    # vector lowvec contains the correct low entries of the matrix.
    #
    if lowvec[cindex] == 0
        return false
    else
        targetindex = lowvec[cindex]
        if psetvec[cindex] == psetvec[targetindex]
            return true
        else
            return false
        end
    end
end

"""
    bool = is_homogeneous(matrix::SparseMatrix, psetvec, cindex)

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

