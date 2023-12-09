export cm_create!, cm_create

"""
    cmatrix, cmatrix_cols = cm_create!(matrix, psetvec)

Compute the connection matrix.

Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`.
"""
function cm_create!(matrix, psetvec)
    #
    # Compute the connection matrix
    #
    numcolumns = size(matrix, 2)
    lowvec = zeros(Int,numcolumns)
    update_low!(matrix, lowvec, startindex=1)

    for j = 1:numcolumns
        jlow = lowvec[j]
        if jlow > 0
            for i = jlow:-1:1
                if matrix[i,j] == 1
                    s = 1
                    found_s = false
                    while (!found_s) & (s <= numcolumns)
                        if ((!(s == j)) & (lowvec[s] == i) & is_homogeneous(matrix, lowvec, psetvec, s))
                            found_s = true
                        else
                            s += 1
                        end
                    end

                    if found_s
                        add_column!(matrix,1,s,j)
                        add_row!(matrix,1,j,s)
                        update_low!(matrix, lowvec, startindex=j)
                    end
                end
            end
        end
    end

    cmatrix_cols = cm_columns(matrix, lowvec, psetvec)
    cmatrix      = matrix[cmatrix_cols, cmatrix_cols]
    return cmatrix, cmatrix_cols
end

"""
    cmatrix, cmatrix_cols, matrixB = cm_create(matrix, psetvec)

Compute the connection matrix.

Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Does not modify the argument `matrix`, but returns
the modified matrix in `matrixB`.
"""
function cm_create(matrix, psetvec)
    #
    # Compute the connection matrix without destroying the input
    #
    matrixB = deepcopy(matrix)
    cmatrix, cmatrix_cols = cm_create!(matrixB, psetvec)
    return cmatrix, cmatrix_cols, matrixB
end

