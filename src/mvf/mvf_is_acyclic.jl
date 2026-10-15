export mvf_is_acyclic

"""
    mvf_is_acyclic(lc::AbstractComplex, mvf::CellSubsets)

Determine whether a multivector field is acyclic.

The function returns a boolean which indicates whether a 
multivector field forms an acyclic partition of the Lefschetz
complex `lc` or not.
"""
function mvf_is_acyclic(lc::AbstractComplex, mvf::CellSubsets)
    #
    # Determine whether a multivector field is acyclic
    #

    # Convert the multivector field to integer form

    if mvf isa Vector{Vector{Int}}
        mvfI = deepcopy(mvf)
    else
        mvfI = convert_cellsubsets(lc, mvf)
    end

    # Append the missing singletons to the mvf

    mvfcells = isempty(mvfI) ? Int[] : reduce(vcat, mvfI)
    critical = setdiff(collect(1:lc.ncells), mvfcells)
    for k in critical
        push!(mvfI, [k])
    end

    # Compute the closures of the multivectors

    mvfIcl = map(t -> lefschetz_closure(lc,t), mvfI)

    # Determine the edges in the directed graph: There is an edge
    # from k to m whenever the closure of multivector k intersects
    # multivector m. These are found via the cell-to-multivector lists.

    cell2mv = cell_to_sets(lc.ncells, mvfI)

    g = DiGraph(length(mvfI))
    for k = 1:length(mvfI)
        for c in mvfIcl[k]
            for m in cell2mv[c]
                if !(k == m)
                    add_edge!(g, k, m)
                end
            end
        end
    end

    # Return the result

    return !is_cyclic(g)
end

