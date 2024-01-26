export homology

"""
    betti = homology(lc::LefschetzComplex; p::Int=2)

Compute the homology of a Lefschetz complex.

The homology is computed over the finite field `GF(p)` and is
returned as a vector `betti` of Betti numbers, where `betti[k]`
is the Betti number in dimension `k-1`.
"""
function homology(lc::LefschetzComplex; p::Int=2)
    #
    # Compute the homology of a Lefschetz complex
    #

    filtration = fill(Int(1),lc.ncells)
    phs, php = persistent_homology(lc,filtration,p=p)
    betti = [length(phs[k]) for k=1:lc.dim+1]
    return betti
end

