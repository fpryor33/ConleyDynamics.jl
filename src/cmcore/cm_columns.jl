export cm_columns

"""
    cm_columns(matrix::SparseMatrix, psetvec::Vector{Int})

Create a vector of column indices for the connection matrix.
"""
function cm_columns(matrix::SparseMatrix, psetvec::Vector{Int})
    #
    # Create a vector of column indices for the connection matrix
    #
    numcolumns = sparse_size(matrix, 2)
    hcols      = homogeneous_columns(matrix, psetvec)
    tcols      = target_columns(matrix, psetvec)
    ccols      = Vector{Int}()

    for j = 1:numcolumns
        if (hcols[j] == false) & (tcols[j] == false)
            push!(ccols,j)
        end
    end

    return ccols
end

