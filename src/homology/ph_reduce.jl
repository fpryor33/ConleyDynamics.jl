export ph_reduce! 

"""
    phsingles, phpairs = ph_reduce!(matrix::SparseMatrix)
    phsingles, phpairs, basis = ph_reduce!(matrix::SparseMatrix;
                                           returnbasis=true)

Apply the persistence reduction algorithm to the matrix.

The function assumes that `matrix` is strictly upper triangular.
It returns the starting columns for infinite length persistence
intervals in `phsingles`, and the birth- and death-columns for
finite length persistence intervals in `phpairs`. If the optional
argument `returnbasis=true` is given, then the function also
returns the computed basis matrix B with `reduced = matrix * B`.
"""
function ph_reduce!(matrix::SparseMatrix; returnbasis::Bool=false)
    #
    # Apply the persistence reduction algorithm to the matrix
    #

    if !sparse_is_sut(matrix)
        error("The matrix has to be strictly upper triangular!")
        return
    end

    # Extract the zero and one elements

    tzero = matrix.zero
    tone  = matrix.one

    # Create the identity matrix for basis computation

    numcolumns = sparse_size(matrix, 2)
    if returnbasis
        basis = sparse_identity(numcolumns, tone)
    end

    # Initialize the main computation

    lowtocolumn = fill(Int(0),numcolumns)
    partofinterval = fill(false,numcolumns)

    for j=1:numcolumns
        keepgoing = true
        while (length(matrix.columns[j]) > 0) & keepgoing
            columnlow = sparse_low(matrix,j)
            if lowtocolumn[columnlow] == 0
                keepgoing = false
                lowtocolumn[columnlow] = j
                partofinterval[columnlow] = true
                partofinterval[j] = true
            else
                s = lowtocolumn[columnlow]
                gamma = matrix[columnlow,j] / matrix[columnlow,s]
                sparse_add_column!(matrix,j,s,-gamma)
                if returnbasis
                    sparse_add_column!(basis,j,s,-gamma)
                end
            end
        end
    end

    # Prepare the return arrays

    phpairs = Vector{Tuple{Int,Int}}()
    phsingles = Vector{Int}()

    for k=1:numcolumns
        if lowtocolumn[k] > 0
            push!(phpairs,(k,lowtocolumn[k]))
        end
        if !partofinterval[k]
            push!(phsingles,k)
        end
    end

    # Return the results

    if returnbasis
        return phsingles, phpairs, basis
    else
        return phsingles, phpairs
    end
end

