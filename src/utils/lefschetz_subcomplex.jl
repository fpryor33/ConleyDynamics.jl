export lefschetz_subcomplex

"""
    lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{Int})

Extract a closed subcomplex from a Lefschetz complex. The subcomplex is
the closure of the collection of cells given in `subcomp`.
"""
function lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{Int})
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    # Determine the closed subcomplex

    clsubcomp = lefschetz_closure(lc, subcomp)

    # Determine the new struct fields

    sub_ncells     = length(clsubcomp)
    sub_labels     = lc.labels[clsubcomp]
    sub_dimensions = lc.dimensions[clsubcomp]
    sub_dim        = maximum(sub_dimensions)
    sub_indices = Dict{String,Int}(sub_labels[j] => j
                                   for j=1:length(sub_labels))

    # Construct the boundary matrix as a sparse matrix

    if typeof(lc.boundary)==Matrix
        bnd = sparse_from_full(lc.boundary)
    else
        bnd = lc.boundary
    end

    sub_boundary = sparse_minor(bnd, clsubcomp, clsubcomp)

    # Create the new Lefschetz complex and return it

    sub_lc = LefschetzComplex(sub_ncells, sub_dim, sub_boundary,
                              sub_labels, sub_indices, sub_dimensions)
    return sub_lc
end

"""
    lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{String})

Extract a closed subcomplex from a Lefschetz complex. The subcomplex is
the closure of the collection of cells given in `subcomp`.
"""
function lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{String})
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcsub = lefschetz_subcomplex(lc, subcompI)
    return lcsub
end

