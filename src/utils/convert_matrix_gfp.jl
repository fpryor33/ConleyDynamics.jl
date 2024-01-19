export convert_matrix_gfp

"""
    gfpmatrix = convert_matrix_gfp(matrix, p)

Convert an integer matrix to a finite field matrix over GF(p).
"""
function convert_matrix_gfp(matrix, p)
    #
    # Convert an integer matrix to a finite field matrix over GF(p).
    #
    m,n = size(matrix)
    FF = GF(p)
    MM = matrix_space(FF, m, n)
    gfpmatrix = MM(matrix)
    return gfpmatrix
end

