export cm_reduce! 

"""
    cm_reduce!(matrix::SparseMatrix, psetvec::Vector{Int};
               [returnbasis::Bool])

Compute the connection matrix.

Assumes that `matrix` is upper triangular and filtered according
to `psetvec`. Modifies the argument `matrix`.

# Return values:
* `cmatrix`: Connection matrix
* `cmatrix_cols`: Columns of the connection matrix in the boundary
* `basisvecs` (optional): If the argument `returnbasis=true` is given,
  this returns information about the computed basis. The k-th entry
  of `basisvecs` is a vector containing the columns making up the
  k-th basis vector, which corresponds to column `cmatrix_cols[k]`.
"""
function cm_reduce!(matrix::SparseMatrix, psetvec::Vector{Int};
                    returnbasis::Bool=false)
    #
    # Compute the connection matrix
    #

    # Extract the zero and one elements

    tchar = matrix.char
    tzero = matrix.zero
    tone  = matrix.one

    # Create the identity matrix for basis computation

    numcolumns = sparse_size(matrix, 2)
    if returnbasis
        basis = sparse_identity(numcolumns, tone, p=tchar)
    end

    # Initialize the main computation

    for j = 1:numcolumns
        jlow = sparse_low(matrix,j)
        if jlow > 0
            for i = jlow:-1:1
                if !(matrix[i,j] == tzero)
                    s = 1
                    found_s = false
                    while (!found_s) & (s <= numcolumns)
                        if ((!(s == j)) & (sparse_low(matrix,s) == i) & is_homogeneous(matrix,psetvec,s))
                            found_s = true
                        else
                            s += 1
                        end
                    end

                    if found_s
                        gamma1 = matrix[i,j]
                        gamma2 = matrix[i,s]
                        sparse_add_column!(matrix,j,s,-gamma1,gamma2)
                        if returnbasis
                            sparse_add_column!(basis,j,s,-gamma1,gamma2)
                        end
                        sparse_add_row!(matrix,s,j,gamma1,gamma2)
                    end
                end
            end
        end
    end

    cmatrix_cols = cm_columns(matrix, psetvec)
    cmatrix      = sparse_minor(matrix, cmatrix_cols, cmatrix_cols)
    
    if returnbasis

        # Create a vector of basis vectors

        basisvecs = Vector{Vector{Int}}()

        for k=1:length(cmatrix_cols)
            bvec = sort(sparse_get_nz_column(basis,cmatrix_cols[k]),rev=true)
            push!(basisvecs,bvec)
        end

        # Return the results

        return cmatrix, cmatrix_cols, basisvecs
    else
        return cmatrix, cmatrix_cols
    end
end

