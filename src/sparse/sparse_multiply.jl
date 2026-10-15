export sparse_multiply

"""
    sparse_multiply(A::SparseMatrix, B::SparseMatrix)

Multiply two sparse matrices.

Exceptions are raised if the matrix product is not defined 
or the entry types do not match.
"""
function sparse_multiply(A::SparseMatrix, B::SparseMatrix)
    #
    # Multiply two sparse matrices.
    #

    # Check whether the product exists

    if !(A.ncol == B.nrow)
        error("Matrix product undefined!")
    end

    if !(typeof(A.zero) == typeof(B.zero))
        error("The sparse matrices need to have the same type!")
    end

    if !(A.char == B.char)
        error("The sparse matrices have to be over the same field!")
    end

    # Column m of the product is the linear combination of the columns
    # of A with the coefficients in column m of B. This only touches
    # the nonzero entries which actually contribute to the product.

    p       = A.char
    tzero   = A.zero
    T       = typeof(tzero)
    entries = [Vector{T}() for _ in 1:B.ncol]
    columns = [Vector{Int}() for _ in 1:B.ncol]

    # The columns are split into chunks, and every chunk uses its own
    # dense accumulator together with the list of rows it touched

    nchunks = min(Threads.nthreads(), B.ncol)
    chunks  = [m:nchunks:B.ncol for m in 1:nchunks]

    Threads.@threads for chunk in chunks
        acc     = fill(tzero, A.nrow)
        touched = fill(false, A.nrow)
        rlist   = Vector{Int}()
        for m in chunk
            for (j, bval) in zip(B.columns[m], B.entries[m])
                for (i, aval) in zip(A.columns[j], A.entries[j])
                    acc[i] = scalar_add(acc[i], scalar_multiply(aval, bval, p), p)
                    if !touched[i]
                        touched[i] = true
                        push!(rlist, i)
                    end
                end
            end
            sort!(rlist)
            for i in rlist
                if !(acc[i] == tzero)
                    push!(columns[m], i)
                    push!(entries[m], acc[i])
                end
                acc[i]     = tzero
                touched[i] = false
            end
            empty!(rlist)
        end
    end

    # Create the row representation, which is automatically sorted

    rows = [Vector{Int}() for _ in 1:A.nrow]
    for m in 1:B.ncol
        for i in columns[m]
            push!(rows[i], m)
        end
    end

    # Construct and return the sparse product matrix

    sm = SparseMatrix{T}(A.nrow, B.ncol, p, tzero, A.one, entries, columns, rows)
    return sm
end

"""
    sparse_multiply(alpha::Int, A::SparseMatrix)

Multiply a scalar with a sparse matrix.

This method is for finite fields. An exception is raised if
the sparse matrix is not defined over a finite field.
"""
function sparse_multiply(alpha::Int, A::SparseMatrix)
    #
    # Multiply a sparse matrix by a scalar
    #
    p = A.char
    @assert p > 0 "Wrong scalar type!"

    if alpha == 0
        B = sparse_zero(A.nrow, A.ncol, p=p)
    else
        B = deepcopy(A)
    end
    for k in 1:A.ncol
        for m in 1:length(A.columns[k])
            val = A.entries[k][m]
            B.entries[k][m] = scalar_multiply(val, alpha, p)
        end
    end
    return B
end

"""
    sparse_multiply(alpha::Rational{Int}, A::SparseMatrix)

Multiply a scalar with a sparse matrix.

This method is for the rational numbers. An exception is raised
if the sparse matrix is not defined over field of rationals.
"""
function sparse_multiply(alpha::Rational{Int}, A::SparseMatrix)
    #
    # Multiply a sparse matrix by a scalar
    #
    p = A.char
    @assert p == 0 "Wrong scalar type!"

    if alpha == 0
        B = sparse_zero(A.nrow, A.ncol, p=p)
    else
        B = deepcopy(A)
    end
    for k in 1:A.ncol
        for m in 1:length(A.columns[k])
            B.entries[k][m] = alpha * A.entries[k][m]
        end
    end
    return B
end

"""
    Base.:*(A::SparseMatrix, B::SparseMatrix)

Multiply two sparse matrices.

Exceptions are raised if the matrix product is not defined 
or the entry types do not match.
"""
Base.:*(A::SparseMatrix,B::SparseMatrix) = sparse_multiply(A,B)

"""
    Base.:*(alpha::Int, A::SparseMatrix)

Multiply a scalar with a sparse matrix.

This method is for finite fields. An exception is raised if
the sparse matrix is not defined over a finite field.
"""
Base.:*(alpha::Int, A::SparseMatrix) = sparse_multiply(alpha,A)

"""
    Base.:*(alpha::Rational{Int}, A::SparseMatrix)

Multiply a scalar with a sparse matrix.

This method is for the rational numbers. An exception is raised
if the sparse matrix is not defined over field of rationals.
"""
Base.:*(alpha::Rational{Int}, A::SparseMatrix) = sparse_multiply(alpha,A)

