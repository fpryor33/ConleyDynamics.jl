export conley_index

"""
    conley_index(lc::LefschetzComplex, subcomp::Vector{String};
                 p::Int)

Determine the Conley index of a Lefschetz complex subset.

The function raises an error if the subset `subcomp` is not
locally closed. The optional parameter `p` specifies the
field characteristic for the homology computation.
"""
function conley_index(lc::LefschetzComplex, subcomp::Vector{String};
                      p::Int=-1)
    #
    # Determine the Conley index of a Lefschetz complex subset
    #

    # Determine the closure-mouth-pair

    clo, mo = lefschetz_clomo_pair(lc, subcomp)

    # Do we have a locally closed set?

    if !lefschetz_is_closed(lc, mo)
        error("This is not a locally closed set!")
    end

    # Compute the subcomplex for the closure

    lcclo = lefschetz_closed_subcomplex(lc, clo)

    # Compute and return the Conley index

    bettisub = relative_homology(lcclo, mo, p=p)
    betti = fill(Int(0),lc.dim+1)
    for k=1:length(bettisub)
        betti[k] = bettisub[k]
    end

    # Return the result

    return betti
end

"""
    conley_index(lc::LefschetzComplex, subcomp::Vector{Int};
                 p::Int)

Determine the Conley index of a Lefschetz complex subset.

The function raises an error if the subset `subcomp` is not
locally closed. The optional parameter `p` specifies the
field characteristic for the homology computation.
"""
function conley_index(lc::LefschetzComplex, subcomp::Vector{Int};
                      p::Int=-1)
    #
    # Determine the Conley index of a Lefschetz complex subset
    #

    subcompS = lc.labels[subcomp]
    return conley_index(lc, subcompS, p=p)
end

