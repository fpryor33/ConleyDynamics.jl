export morse_sets

"""
    morse_sets(lc::LefschetzComplex, mvf::CellSubsets; poset::Bool=false)

Find the nontrivial Morse sets of a multivector field on a Lefschetz complex.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The function returns the nontrivial Morse sets as
a `Vector{Vector{Int}}`. If the optional argument `poset=true1` is added,
then the function returns both the Morse sets and the adjacency matrix of
the Hasse diagram of the underlying poset.
"""
function morse_sets(lc::LefschetzComplex, mvf::CellSubsets; poset::Bool=false)
    #
    # Find the nontrivial Morse sets of a multivector field
    #

    # Convert the multivector field to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = mvf
        return_as_labels = false
    else
        mvfI = convert_cellsubsets(lc, mvf)
        return_as_labels = true
    end
    
    # Create the digraph based on the boundary matrix

    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(lc.boundary)
    edgelist = Edge.([(c[k],r[k]) for k in 1:length(c)])
    dg = SimpleDiGraph(edgelist)

    # Add cliques for the multivectors

    lenmvf = length(mvfI)
    for k = 1:lenmvf
        mv  = mvfI[k]
        lmv = length(mv)
        for m = 1:lmv-1
            for n = m+1:lmv
                add_edge!(dg,mv[m],mv[n])
                add_edge!(dg,mv[n],mv[m])
            end
        end
    end

    # Apply Tarjan's algorithm
    
    scc = strongly_connected_components_tarjan(dg)

    for k=1:length(scc)
        sort!(scc[k])
    end

    # Find all Morse sets with nontrivial homology

    isnontrivial = fill(Int(0), length(scc))

    Threads.@threads for k = 1:length(scc)
        isnontrivial[k] = sum(conley_index(lc, scc[k]))
    end

    ntindices = findall(x -> x>0, isnontrivial)
    morsesets = scc[ntindices]

    # Determine the Hasse diagram of the poset

    if poset
        nms = length(morsesets)
        psadj = fill(false,nms,nms)
        for k1 in 1:nms
            vk1 = morsesets[k1][1]
            for k2 in 1:nms
                vk2 = morsesets[k2][1]
                psadj[k1,k2] = has_path(dg, vk1, vk2)
            end
        end
        psdg = SimpleDiGraph(psadj)
        pshasse = transitivereduction(psdg)

        hasse = fill(false,nms,nms)
        for k1 in 1:nms
            for k2 in 1:nms
                hasse[k2,k1] = has_edge(pshasse, k1, k2)
            end
        end
    end

    # Return the results

    if return_as_labels
        morsesets = convert_cellsubsets(lc, morsesets)
    end

    if poset
        return morsesets, hasse
    else
        return morsesets
    end
end

