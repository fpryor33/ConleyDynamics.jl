export lefschetz_closure

"""
    c = lefschetz_closure(lc::LefschetzComplex, subcomp::Vector{Int})

Compute the closure of a Lefschetz complex subset.
"""
function lefschetz_closure(lc::LefschetzComplex, subcomp::Vector{Int})
    #
    # Compute the closure of a Lefschetz complex subset
    #

    # Extract the boundary matrix as a sparse matrix

    if typeof(lc.boundary)==Matrix
        bnd = sparse_from_full(lc.boundary)
    else
        bnd = lc.boundary
    end

    # Create list of cells to consider

    lcsubset = deepcopy(subcomp)
    lcclosure = Vector{Int}()

    while length(lcsubset)>0
        cell = popfirst!(lcsubset)
        push!(lcclosure, cell)
        if length(bnd.columns[cell])>0
            append!(lcsubset, bnd.columns[cell])
        end
    end

    # Remove duplicates, sort, and return

    return sort(unique(lcclosure))
end

"""
    c = lefschetz_closure(lc::LefschetzComplex, subcomp::Vector{String})

Compute the closure of a Lefschetz complex subset.
"""
function lefschetz_closure(lc::LefschetzComplex, subcomp::Vector{String})
    #
    # Compute the closure of a Lefschetz complex subset
    #

    subcompI = Vector{Int}()
    for ks in subcomp
        push!(subcompI,lc.indices[ks])
    end

    lcclosure = lefschetz_closure(lc, subcompI)
    return lcclosure
end

