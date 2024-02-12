export homology

"""
    homology(lc::LefschetzComplex; [p::Int])

Compute the homology of a Lefschetz complex.

The homology is computed over the rationals (for `p=0`) or the
finite field `GF(p)` (for prime `p`) and is returned as a vector
`betti` of Betti numbers, where `betti[k]` is the Betti number
in dimension `k-1`. If the Lefschetz complex boundary matrix
already has been specialized to a field, the optional
argument `p` can be omitted.
"""
function homology(lc::LefschetzComplex; p::Int=-1)
    #
    # Compute the homology of a Lefschetz complex
    #

    filtration = fill(Int(1),lc.ncells)
    phs, php = persistent_homology(lc,filtration,p=p)
    betti = [length(phs[k]) for k=1:lc.dim+1]
    return betti
end

