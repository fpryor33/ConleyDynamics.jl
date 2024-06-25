export lefschetz_subcomplex

"""
    lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{Int})

Extract a subcomplex from a Lefschetz complex. The subcomplex has to be
locally closed, and is given by the collection of cells in `subcomp`.
"""
function lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{Int})
    #
    # Extract a subcomplex of a Lefschetz complex
    #

    # Do we have a locally closed set?

    if !lefschetz_is_locally_closed(lc, subcomp)
        error("This is not a locally closed set!")
    end

    # Determine the new struct fields

    subcompsorted  = sort(subcomp)
    sub_ncells     = length(subcompsorted)
    sub_labels     = lc.labels[subcompsorted]
    sub_dimensions = lc.dimensions[subcompsorted]
    sub_dim        = maximum(sub_dimensions)
    sub_indices = Dict{String,Int}(sub_labels[j] => j
                                   for j=1:length(sub_labels))

    # Construct the new boundary matrix

    bnd = lc.boundary
    sub_boundary = sparse_minor(bnd, subcompsorted, subcompsorted)

    # Create the new Lefschetz complex and return it

    sub_lc = LefschetzComplex(sub_ncells, sub_dim, sub_boundary,
                              sub_labels, sub_indices, sub_dimensions)
    return sub_lc
end

"""
    lefschetz_subcomplex(lc::LefschetzComplex, subcomp::Vector{String})

Extract a subcomplex from a Lefschetz complex. The subcomplex has to be
locally closed, and is given by the collection of cells in `subcomp`.
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

