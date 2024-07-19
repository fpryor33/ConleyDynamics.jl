export morse_sets

"""
    morse_sets(lc::LefschetzComplex, mvf::CellSubsets)

Find the nontrivial Morse sets of a multivector field on a Lefschetz complex.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The function returns the nontrivial Morse sets as
a `Vector{Vector{Int}}`.
"""
function morse_sets(lc::LefschetzComplex, mvf::CellSubsets)
    #
    # Find the nontrivial Morse sets of a multivector field
    #

    # Convert the multivector field to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = mvf
    else
        mvfI = convert_cellsubsets(lc, mvf)
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

    # Return the results

    return morsesets
end

