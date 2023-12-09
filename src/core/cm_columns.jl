export cm_columns

"""
    ccols = cm_columns(matrix, lowvec, psetvec)

Create a vector of column indices for the connection matrix.
"""
function cm_columns(matrix, lowvec, psetvec)
    #
    # Create a vector of column indices for the connection matrix
    #
    numcolumns = size(matrix, 2)
    hcols      = homogeneous_columns(matrix, lowvec, psetvec)
    tcols      = target_columns(matrix, lowvec, psetvec)
    ccols      = Vector{Int}()

    for j = 1:numcolumns
        if (hcols[j] == false) & (tcols[j] == false)
            push!(ccols,j)
        end
    end

    return ccols
end

