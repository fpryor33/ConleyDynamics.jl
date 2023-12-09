export target_columns
"""
    tcols = target_columns(matrix, lowvec, psetvec)

Determine which columns of `matrix` are target columns.
"""
function target_columns(matrix, lowvec, psetvec)
    #
    # Determine which columns are target columns
    #
    numcolumns = size(matrix, 2)
    tcols      = Vector(zeros(Bool,numcolumns))

    for j = 1:numcolumns
        if is_homogeneous(matrix, lowvec, psetvec, j)
            tcols[lowvec[j]] = true
        end
    end

    return tcols
end

