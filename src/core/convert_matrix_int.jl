export convert_matrix_int

"""
    intmatrix = convert_matrix_int(matrix)

Convert a matrix over the finite field GF(p) to an integer matrix.
"""
function convert_matrix_int(matrix)
    #
    # Convert a matrix with finite field entries to an integer matrix
    #
    m,n = size(matrix)
    intmatrix = [Int(lift(matrix[i,j])) for i = 1:m, j = 1:n]
    return intmatrix
end

