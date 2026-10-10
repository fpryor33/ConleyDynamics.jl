export lefschetz_subcomplex

"""
    lefschetz_subcomplex(lc::AbstractComplex, subcomp::Vector{Int};
                         [validate::Bool=true])

Extract a subcomplex from a Lefschetz complex. The subcomplex has to be
locally closed, and is given by the collection of cells in `subcomp`.

If `lc` is a `EuclideanComplex`, the returned subcomplex is also an
`EuclideanComplex` with the coordinates restricted to the subcomplex cells.

By default, the function raises an error if `subcomp` is not locally
closed. This check can be skipped by passing the optional argument
`validate=false`, which should only be done if it is already known
that the set is locally closed. The computational cost of the function
is proportional to the size of `subcomp`, not to the size of `lc`.
"""
function lefschetz_subcomplex(lc::AbstractComplex, subcomp::Vector{Int};
                              validate::Bool=true)
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    # Do we have a locally closed set?

    if validate && !lefschetz_is_locally_closed(lc, subcomp)
        error("This is not a locally closed set!")
    end

    # Determine the new struct fields

    subcompsorted  = sort(subcomp)
    sub_labels     = lc.labels[subcompsorted]
    sub_dimensions = lc.dimensions[subcompsorted]

    # Construct the new boundary matrix

    sub_boundary = lefschetz_subboundary(lc, subcompsorted)

    # Create the new complex and return it

    if lc isa EuclideanComplex
        sub_coords = lc.coords[subcompsorted]
        return EuclideanComplex(sub_labels, sub_dimensions, sub_boundary,
                                sub_coords; validate=false)
    else
        return LefschetzComplex(sub_labels, sub_dimensions, sub_boundary; validate=false)
    end
end

"""
    lefschetz_subcomplex(lc::AbstractComplex, subcomp::Vector{String};
                         [validate::Bool=true])

Extract a subcomplex from a Lefschetz complex. The subcomplex has to be
locally closed, and is given by the collection of cells in `subcomp`.

If `lc` is a `EuclideanComplex`, the returned subcomplex is also an
`EuclideanComplex` with the coordinates restricted to the subcomplex cells.

By default, the function raises an error if `subcomp` is not locally
closed. This check can be skipped by passing the optional argument
`validate=false`, which should only be done if it is already known
that the set is locally closed.
"""
function lefschetz_subcomplex(lc::AbstractComplex, subcomp::Vector{String};
                              validate::Bool=true)
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcsub = lefschetz_subcomplex(lc, subcompI, validate=validate)
    return lcsub
end

#
# lefschetz_subboundary(lc::AbstractComplex, cells::Vector{Int})
#
# Restrict the boundary matrix of a Lefschetz complex to a set of cells.
#
# This internal helper function returns the minor of `lc.boundary` whose
# rows and columns are given by `cells`, which has to be sorted and cannot
# contain repetitions. In contrast to `sparse_minor`, it only accesses the
# columns of the boundary matrix which correspond to `cells`, and its
# computational cost is therefore independent of the size of `lc`.
#
function lefschetz_subboundary(lc::AbstractComplex, cells::Vector{Int})
    #
    # Restrict the boundary matrix to a sorted list of cells
    #

    bnd = lc.boundary
    T   = typeof(bnd.zero)
    n   = length(cells)

    entries = [Vector{T}() for _ in 1:n]
    columns = [Vector{Int}() for _ in 1:n]
    rows    = [Vector{Int}() for _ in 1:n]

    # The map from cells to the new indices is increasing, so
    # the resulting column and row index vectors are sorted

    for j in 1:n
        bcol = bnd.columns[cells[j]]
        bent = bnd.entries[cells[j]]
        for m in eachindex(bcol)
            i = searchsortedfirst(cells, bcol[m])
            if (i <= n) && (cells[i] == bcol[m])
                push!(columns[j], i)
                push!(entries[j], bent[m])
                push!(rows[i], j)
            end
        end
    end

    return SparseMatrix{T}(n, n, bnd.char, bnd.zero, bnd.one,
                           entries, columns, rows)
end
