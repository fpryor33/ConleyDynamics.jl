export admissible_order

"""
    admiorder, admibnd, sccnumber, scc = admissible_order(bndmatrix, mvf)

Find an admissible order based on the boundary matrix and the multivector field.

The vector `sccnumber` contains the strongly connected component for each column,
and the matrix `admibnd` is the reordered boundary matrix corresponding to the
order given by `admiorder`.
"""
function admissible_order(bndmatrix::Matrix{Int}, mvf::Vector{Vector{Int}})
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

    # Create the reordered boundary matrix

    admibnd = bndmatrix[admiorder,admiorder]

    # Return the results

    return admiorder, admibnd, sccnumber, scc
end

"""
    admiorder, admibnd, sccnumber, scc = admissible_order(bndmatrix, mvf)

Find an admissible order based on the boundary matrix and the multivector field.

The vector `sccnumber` contains the strongly connected component for each column,
and the matrix `admibnd` is the reordered boundary matrix corresponding to the
order given by `admiorder`.
"""
function admissible_order(bndmatrix::SparseMatrix{Int}, mvf::Vector{Vector{Int}})
    #
    # Find an admissible order based on the boundary matrix and the multivector
    # field
    #
    
    # Create the digraph based on the boundary matrix

    nr, nc, tzero, tone, r, c, vals = lists_from_sparse(bndmatrix)
    edgelist = Edge.([(c[k],r[k]) for k in 1:length(c)])
    dg = SimpleDiGraph(edgelist)

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

    # Create the reordered boundary matrix

    admibnd = sparse_permute(bndmatrix, admiorder, admiorder)

    # Return the results

    return admiorder, admibnd, sccnumber, scc
end

