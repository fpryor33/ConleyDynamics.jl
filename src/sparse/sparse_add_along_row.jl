export sparse_add_along_row!

"""
    sparse_add_along_row!(matrix::SparseMatrix, ri::Int, ci::Int)

Remove the sparse matrix entry at location `(ri,ci)`.
"""
function sparse_add_along_row!(matrix::SparseMatrix, ri::Int, ci1::Int, ci2::Int, c)
    #
    # Replace matrix[ri,ci1] by matrix[ri,ci1] + c * matrix[ri,ci2]
    #

    # Find the location of (ri,ci2) in column ci2

    index_in_col2 = findfirst(x -> x==ri, matrix.columns[ci2])

    # If matrix[ri,ci2] == 0 do nothing

    if index_in_col2 == nothing
        return
    end

    # The element is nonzero, so we have to do something..
    # First, determine the new value of matrix[ri,ci1].

    index_in_col1 = findfirst(x -> x==ri, matrix.columns[ci1])
    new_value1 = c * matrix.entries[ci2][index_in_col2]

    if !(index_in_col1 == nothing)
        new_value1 += matrix.entries[ci1][index_in_col1]
    end

    # Second, incorporate the updated value into the matrix

    if new_value1 == matrix.zero
        # Entry became 0, has to be removed, if present
        if !(index_in_col1 == nothing)
            sparse_remove!(matrix, ri, ci1)
        end
    elseif !(index_in_col1 == nothing)
        # Entry present, just needs to be update
        matrix.entries[ci1][index_in_col1] = new_value1
    else
        # Insert ci1 into row ri
        insert_index_row = findfirst(x -> x>ci1, matrix.rows[ri])
        if insert_index_row == nothing
            push!(matrix.rows[ri], ci1)
        else
            insert!(matrix.rows[ri], insert_index_row, ci1)
        end

        # Insert ri into column ci1 and add the new entry
        insert_index_col = findfirst(x -> x>ri, matrix.columns[ci1])
        if insert_index_col == nothing
            push!(matrix.columns[ci1], ri)
            push!(matrix.entries[ci1], new_value1)
        else
            insert!(matrix.columns[ci1], insert_index_col, ri)
            insert!(matrix.entries[ci1], insert_index_col, new_value1)
        end
    end

    return
end

