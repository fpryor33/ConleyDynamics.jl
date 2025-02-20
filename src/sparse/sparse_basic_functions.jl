export sparse_size, sparse_low, sparse_identity
export sparse_fullness, sparse_sparsity
export sparse_is_zero
export sparse_show

"""
    sparse_size(matrix::SparseMatrix, dim::Int)

Number of rows (`dim=1`) or columns (`dim=2`) of a sparse matrix.
"""
function sparse_size(matrix::SparseMatrix, dim::Int)
    #
    # Return the number of columns or rows
    #

    if dim == 1
        return matrix.nrow
    elseif dim == 2
        return matrix.ncol
    else
        error("dim has to be 1 or 2!")
    end
end

"""
    sparse_low(matrix::SparseMatrix, col::Int)

Row index of the lowest nonzero matrix entry in column `col`.
"""
function sparse_low(matrix::SparseMatrix, col::Int)
    #
    # Return the location of the lowest nonzero entry in 
    # column col of the matrix
    #

    if length(matrix.columns[col]) == 0
        return 0
    else
        return last(matrix.columns[col])
    end
end

"""
    sparse_identity(n::Int; p::Int=0)

Create a sparse identity matrix with n rows and columns.

The optional argument `p` specifies the field characteristic.
If `p=0` then the sparse matrix is over the rationals, while
if `p>0` is a prime, then the matrix is an integer matrix 
whose entries are interpreted in `GF(p)`.
"""
function sparse_identity(n::Int; p::Int=0)
    #
    # Create a sparse identity matrix with n rows and columns,
    # and with diagonal entry tone.
    #

    # Initialize the variables and lists

    if p == 0
        tzero = Rational{Int}(0)
        tone  = Rational{Int}(1)
    elseif p > 0
        tzero = Int(0)
        tone  = Int(1)
    else
        error("The characteristic cannot be negative!")
    end

    r = Vector{Int}(1:n)
    vals = fill(tone,n)

    # Create the sparse identity and return it

    sm = sparse_from_lists(n, n, p, tzero, tone, r, r, vals)
    return sm
end

"""
    sparse_fullness(sm::SparseMatrix)

Display the fullness of the sparse matrix `sm`.
"""
function sparse_fullness(sm::SparseMatrix)
    #
    # Display the fullness of a sparse matrix
    #

    nz = 0
    for k=1:sm.ncol
        nz += length(sm.columns[k])
    end
    fullness = Float64(nz) / (Float64(sm.ncol) * Float64(sm.nrow))
    
    return fullness
end

"""
    sparse_sparsity(sm::SparseMatrix)

Display the sparsity of the sparse matrix `sm`.
"""
function sparse_sparsity(sm::SparseMatrix)
    #
    # Display the sparsity of a sparse matrix
    #

    return 1.0 - sparse_fullness(sm)
end

"""
    sparse_is_zero(sm::SparseMatrix)

Test whether the sparse matrix `sm` is the zero matrix.
"""
function sparse_is_zero(sm::SparseMatrix)
    #
    # Test whether a sparse matrix is zero
    #

    nz = 0
    for k=1:sm.ncol
        nz += length(sm.columns[k])
    end

    if nz == 0
        return true
    else
        return false
    end
end

"""
    sparse_show(sm::SparseMatrix)

Display the sparse matrix `sm`.
"""
function sparse_show(sm::SparseMatrix)
    #
    # Display a sparse matrix
    #

    for k=1:sm.nrow
        print("[", sm[k,1])
        for m=2:sm.ncol
            print("   ", sm[k,m])
        end
        println("]")
    end
end

# """
#     Base.show(io::IO, sm::SparseMatrix)
# 
# Display the sparse matrix `sm`.
# """
# function Base.show(io::IO, sm::SparseMatrix)
#     #
#     # Display a sparse matrix
#     #
# 
#     for k=1:sm.nrow
#         print(io, "\n[", sm[k,1])
#         for m=2:sm.ncol
#             print(io, "   ", sm[k,m])
#         end
#         print(io, "]")
#     end
# end

