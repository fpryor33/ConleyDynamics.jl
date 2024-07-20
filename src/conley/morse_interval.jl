export morse_interval

"""
    morse_interval(lc::LefschetzComplex, mvf::CellSubsets,
                   ms::CellSubsets)

Find the isolated invariant set for a Morse set interval.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The collection of Morse sets are contained in`ms`.
All of these sets have to be Morse sets in the sense of being strongly
connected components of the flow graph. (For runtime issues, this will
not be checked!) In other words, the sets in `ms` have to be determined
using the function `morse_sets`!

The function returns the smallest isolated invariant set which contains
the Morse sets and their connections as a `Vector{Int}`.
"""
function morse_interval(lc::LefschetzComplex, mvf::CellSubsets,
                        ms::CellSubsets)
    #
    # Find the isolated invariant set for a Morse set interval
    #

    # Convert the multivector field and the Morse sets to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = mvf
    else
        mvfI = convert_cellsubsets(lc, mvf)
    end

    if ms isa Vector{Vector{Int}}
        msI = ms
    else
        msI = convert_cellsubsets(lc, ms)
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

    # Add cliques for the Morse sets
    #
    #  THIS IS NOT NECESSARY IF EACH MORSE SET
    #  IS A STRONGLY CONNECTED COMPONENT!!
    #
    #  That's why it is commented out..
    #
    # lenms = length(msI)
    # for k = 1:lenms
    #     cs = msI[k]
    #     lcs = length(cs)
    #     for m = 1:lcs-1
    #         for n = m+1:lcs
    #             add_edge!(dg,cs[m],cs[n])
    #             add_edge!(dg,cs[n],cs[m])
    #         end
    #     end
    # end

    # Connect the Morse sets to each other

    for k=2:length(msI)
        add_edge!(dg,msI[1][1],msI[k][1])
        add_edge!(dg,msI[k][1],msI[1][1])
    end

    # Apply Tarjan's algorithm
    
    scc = strongly_connected_components_tarjan(dg)

    # Find the strongly connected component which contains
    # the Morse sets in ms

    invhullindex = first(findall(x -> msI[1][1] in x, scc))
    invhull = sort(scc[invhullindex])

    # Return the results

    return invhull
end

