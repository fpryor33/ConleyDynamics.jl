export admissible_order

"""
    admiorder, sccnumber, scc = admissible_order(bndmatrix, mvf)

Find an admissible order based on the boundary matrix and the multivector field.

The vector sccnumber contains the strongly connected component for each column.
"""
function admissible_order(bndmatrix, mvf)
    #
    # Find an admissible order based on the boundary matrix and the multivector
    # field
    #
    
    # Create the digraph based on the boundary matrix

    dg = SimpleDiGraph(bndmatrix')

    # Add cliques for the multivectors

    lenmvf = length(mvf)
    for k = 1:lenmvf
        mv  = mvf[k]
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

    # Create the admissible order and the vector of scc numbers

    admiorder = Vector{Int}()
    sccnumber = Vector{Int}()

    scccount = 0

    for k=1:length(scc)
        scccount += 1
        append!(admiorder,scc[k])
        append!(sccnumber,scc[k] .* 0 .+ scccount)
    end

    # Return the results

    return admiorder, sccnumber, scc
end

