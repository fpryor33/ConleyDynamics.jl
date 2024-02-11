export sparse_permute

"""
    smp = sparse_permute(sm::SparseMatrix, pr::Vector{Int}, pc::Vector{Int})

Create sparse matrix by permuting the row and column indices.

The vector `pr` describes the row permutation, and `pc` the
column permutation.
"""
function sparse_permute(sm::SparseMatrix, pr::Vector{Int}, pc::Vector{Int})
    #
    # Create sparse matrix by permuting the row and column indices
    #

    # Convert the matrix to list format

    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(sm)

    # Compute the inverse permutations

    ipr = invperm(pr)
    ipc = invperm(pc)

    # Create the permuted index lists

    newr = ipr[r]
    newc = ipc[c]

    # Create and return new sparse matrix

    smp = sparse_from_lists(nr, nc, tchar, tzero, tone, newr, newc, vals)
    return smp
end

