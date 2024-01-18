export sparse_add_column!, sparse_add_row!

"""
    sparse_add_column!(matrix::SparseMatrix, ci1::Int, ci2::Int, c)

Replace `column[ci1]` by `column[ci1] + c * column[ci2]`.
"""
function sparse_add_column!(matrix::SparseMatrix, ci1::Int, ci2::Int, c)
    #
    # Replace column[ci1] by column[ci1] + c * column[ci2]
    #

    # Find all nonzero entries in column[ci2]

    nonzero_rows    = matrix.columns[ci2]
    nonzero_entries = matrix.entries[ci2]

    # Loop though the entries and perform the column operation

    for k in 1:length(nonzero_rows)
        row_index = nonzero_rows[k]
        new_value = c * nonzero_entries[k]
        new_value += sparse_get_entry(matrix, row_index, ci1)
        sparse_set_entry!(matrix, row_index, ci1, new_value)
    end
end

"""
    sparse_add_row!(matrix::SparseMatrix, ri1::Int, ri2::Int, c)

Replace `row[ri1]` by `row[ri1] + c * row[ri2]`.
"""
function sparse_add_row!(matrix::SparseMatrix, ri1::Int, ri2::Int, c)
    #
    # Replace row[ri1] by row[ri1] + c * row[ri2]
    #

    # Find all nonzero entries in row[ri2]

    nonzero_columns = matrix.rows[ri2]

    # Loop though the entries and perform the column operation

    for kc in nonzero_columns
        new_value = c * sparse_get_entry(matrix, ri2, kc)
        new_value += sparse_get_entry(matrix, ri1, kc)
        sparse_set_entry!(matrix, ri1, kc, new_value)
    end
end

