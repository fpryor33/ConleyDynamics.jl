export create_simplicial_complex

"""
    lc = create_simplicial_complex(labels::Vector{String},
                                   simplices::Vector{Vector{Int}})

Initialize a Lefschetz complex from a simplicial complex.

The vector `labels` contains a label for every vertex, while
`simplices` contains all the highest-dimensional simplices necessary
to define the simplicial complex.
"""
function create_simplicial_complex(labels::Vector{String},
                                   simplices::Vector{Vector{Int}})
    #
    # Create a Lefschetz complex struct for a simplicial complex.
    #
    
    sclabels = deepcopy(labels)
    scsimplices = deepcopy(simplices)

    # Find the dimension of the simplicial complex

    scdim = 0
    for k = 1:length(scsimplices)
        scdim = max(scdim, length(scsimplices[k])-1)
    end

    # Create labels for all simplices of dimension at least 1
    # in the simplicial complex and organize them in sets by
    # dimension. In addition, create a dictionary which relates
    # the face labels to the index sequences.
    
    labelsets = Vector{Set{String}}()
    for k=1:scdim
        push!(labelsets,Set{String}())
    end

    labelvertexdict = Dict{String,Vector{Int}}()

    for k=1:length(scsimplices)
        csimp = scsimplices[k]
        cdim = length(csimp) - 1
        for m=2:cdim+1
            for faceindices in combinations(csimp,m)
                facelab = join(sort(sclabels[faceindices]))
                push!(labelsets[m-1],facelab)
                labelvertexdict[facelab] = sort(faceindices)
            end
        end
    end
    
    # Order the simplex labels and create label to index dictionary

    labelsbydim = Vector{Vector{String}}()
    for k=0:scdim
        push!(labelsbydim,Vector{String}())
    end

    append!(labelsbydim[1],sclabels)
    for k=1:scdim
        append!(labelsbydim[k+1],sort(collect(labelsets[k])))
    end

    labelsvec = reduce(vcat,labelsbydim)
    labelindexdict = Dict{String,Int}(labelsvec[j] => j for j=1:length(labelsvec))

    # Create the Poincare polynomials for the simplices

    PP, t = ZZ["t"]
    ppvec = Vector{typeof(t)}()
    for k=0:scdim
        for m=1:length(labelsbydim[k+1])
            push!(ppvec,t^k)
        end
    end

    # Create the boundary map
    
    nsimp  = length(labelsvec)
    nsimp0 = length(sclabels)

    B = zeros(Int,nsimp,nsimp)
    for k=nsimp0+1:nsimp
        csimp = labelvertexdict[labelsvec[k]]
        coeff = 1
        for m=1:length(csimp)
            csimptmp = deepcopy(csimp)
            deleteat!(csimptmp,m)
            bsimp = labelindexdict[join(sclabels[csimptmp])]
            B[bsimp,k] = coeff
            coeff = -coeff
        end
    end
    
    # Create the Lefschetz complex

    lc = LefschetzComplex{typeof(t)}(nsimp,B,labelsvec,labelindexdict,ppvec)

    # Return the Lefschetz complex

    return lc
end

"""
    lc = create_simplicial_complex(labels::Vector{String},
                                   simplices::Vector{Vector{String}})

Initialize a Lefschetz complex from a simplicial complex.

The vector `labels` contains a label for every vertex, while
`simplices` contains all the highest-dimensional simplices necessary
to define the simplicial complex.
"""
function create_simplicial_complex(labels::Vector{String},
                                   simplices::Vector{Vector{String}})
    #
    # Create a Lefschetz complex struct for a simplicial complex.
    #
    # This is an alternative method for simplices represented
    # in String format.
    #
    newsimplices = convert_simplices(simplices, labels)
    lc = create_simplicial_complex(labels, newsimplices)
    return lc
end

