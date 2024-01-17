export sparse_remove!

"""
    sparse_remove!(matrix::SparseMatrix, ri::Int, ci::Int)

Remove the sparse matrix entry at location `(ri,ci)`.
"""
function sparse_remove!(matrix::SparseMatrix, ri::Int, ci::Int)
    #
    # Remove the sparse matrix entry at location (ri,ci).
    #

    # Find the location of (ri,ci) in column ci

    index_in_col = findfirst(x -> x==ri, matrix.columns[ci])

    if !(index_in_col == nothing)
        index_in_row = findfirst(x -> x==ci, matrix.rows[ri])
        deleteat!(matrix.rows[ri],index_in_row)
        deleteat!(matrix.columns[ci],index_in_col)
        deleteat!(matrix.entries[ci],index_in_col)
    end
end

