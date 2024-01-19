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
function cm_create!(matrix, psetvec; returnbasis=false)
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
                if !(matrix[i,j] == 0)
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
                        gamma = matrix[i,j] / matrix[i,s]
                        @views matrix[:,j] = matrix[:,j] .- gamma .* matrix[:,s]
                        if returnbasis
                            @views basis[:,j] = basis[:,j] .- gamma .* basis[:,s]
                        end
                        @views matrix[s,:] = matrix[s,:] .+ gamma .* matrix[j,:]
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

