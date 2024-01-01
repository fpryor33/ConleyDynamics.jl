export cm_create! 

"""
    cmatrix, cmatrix_cols = cm_create!(matrix, psetvec)
    cmatrix, cmatrix_cols, basis = cm_create!(matrix, psetvec;
                                              returnbasis=true)

Compute the connection matrix.

Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`. If the optional
argument `returnbasis=true` is given, then the function also
returns the basis transformation matrix. The columns of `basis`
are the basis vectors which lead to the reduced matrix
representation.
"""
function cm_create!(matrix, psetvec; returnbasis=False)
    #
    # Compute the connection matrix
    #

    # Create the identity matrix for basis computation

    numcolumns = size(matrix, 2)
    if returnbasis
        basis = deepcopy(matrix)
        basis = 0 * basis
        for i=1:numcolumns
            basis[i,i] = 1
        end
    end

    # Initialize the main computation

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
                        if returnbasis
                            add_column!(basis,1,s,j)
                        end
                        add_row!(matrix,1,j,s)
                        update_low!(matrix, lowvec, startindex=j)
                    end
                end
            end
        end
    end

    cmatrix_cols = cm_columns(matrix, lowvec, psetvec)
    cmatrix      = matrix[cmatrix_cols, cmatrix_cols]
    
    if returnbasis
        return cmatrix, cmatrix_cols, basis
    else
        return cmatrix, cmatrix_cols
    end
end

