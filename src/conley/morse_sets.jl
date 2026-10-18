export morse_sets

"""
    morse_sets(lc::AbstractComplex, mvf::CellSubsets; poset::Bool=false,
               poset_algorithm::Symbol=:pairwise)

Find the nontrivial Morse sets of a multivector field on a Lefschetz complex.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The function returns the nontrivial Morse sets as
a `Vector{Vector{Int}}`. If the optional argument `poset=true` is added,
then the function returns both the Morse sets and the adjacency matrix of
the Hasse diagram of the underlying poset.

`poset_algorithm` selects how the Hasse diagram is computed when
`poset=true`:
* `:pairwise` (default) issues one `has_path` query per ordered pair of
  Morse sets (`nms^2` total), each a fresh traversal of the full
  cell-level digraph. Simple, well-tested, but can dominate runtime once
  `nms` grows into the thousands.
* `:condensation` computes the same adjacency matrix via one
  single-source reachability traversal per Morse set (`nms` total,
  parallelized across threads) instead of `nms^2` -- the reachable set
  from each representative is computed once and reused for every
  membership check against the other representatives. Produces
  identical results to `:pairwise`, just asymptotically cheaper for
  large `nms`.
"""
function morse_sets(lc::AbstractComplex, mvf::CellSubsets; poset::Bool=false,
                    poset_algorithm::Symbol=:pairwise)
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

    # Create the cell-to-mv map

    cell2mv = fill(Int(0), lc.ncells)

    for k = 1:length(mvfI)
        for m in mvfI[k]
            if cell2mv[m] == 0
                cell2mv[m] = k
            else
                error(" The multivector field is not a partition!")
            end
        end
    end

    # Create the digraph based on the boundary matrix

    nr, nc, tchar, tzero, tone, r, c, vals = lists_from_sparse(lc.boundary)

    # Add cycles for the multivectors to the lists c and r

    lenmvf = length(mvfI)
    for k = 1:lenmvf
        mv  = mvfI[k]
        lmv = length(mv)
        for m = 1:lmv-1
            push!(c,mv[m])
            push!(r,mv[m+1])
        end
        push!(c,mv[lmv])
        push!(r,mv[1])
    end

    # Create the edge list and the digraph

    edgelist = Edge.([(c[k],r[k]) for k in 1:length(c)])
    dg = SimpleDiGraph(edgelist)

    # Apply Tarjan's algorithm
    
    scc = strongly_connected_components_tarjan(dg)

    for k=1:length(scc)
        sort!(scc[k])
    end

    # Find all nontrivial Morse sets:
    # This means either nontrivial homology or a
    # Morse set containing at least two multivectors

    isnontrivial = fill(Int(0), length(scc))

    # A strongly connected component which contains at least two
    # multivectors is always nontrivial. Otherwise, it is a single
    # multivector, which is locally closed, and its Conley index
    # can be determined from the multivector cells alone

    Threads.@threads for k = 1:length(scc)
        if length(unique(cell2mv[scc[k]])) > 1
            isnontrivial[k] = 1
        else
            isnontrivial[k] = sum(conley_index(lc, scc[k], algorithm="rank",
                                               validate=false))
        end
    end

    # Finally extract the Morse sets

    ntindices = findall(x -> x>0, isnontrivial)
    morsesets = scc[ntindices]

    # Determine the Hasse diagram of the poset

    if poset
        nms = length(morsesets)
        reps = [morsesets[k][1] for k in 1:nms]
        psadj = fill(false,nms,nms)

        if poset_algorithm == :pairwise
            for k1 in 1:nms
                vk1 = reps[k1]
                for k2 in 1:nms
                    vk2 = reps[k2]
                    psadj[k1,k2] = has_path(dg, vk1, vk2)
                end
            end
        elseif poset_algorithm == :condensation
            # One single-source reachability traversal per Morse set
            # (parallelized) instead of nms^2 individual has_path calls --
            # the reachable set from each representative is computed once
            # and reused for every membership check, instead of re-deriving
            # it from scratch for every (k1,k2) pair. Manual BFS with an
            # index-based queue (not popfirst!) to avoid O(n)-per-dequeue
            # overhead on a plain Vector.
            nv_dg = nv(dg)
            Threads.@threads for k1 in 1:nms
                visited = falses(nv_dg)
                queue = Vector{Int}(undef, nv_dg)
                queue[1] = reps[k1]
                visited[reps[k1]] = true
                qhead = 1
                qtail = 1
                while qhead <= qtail
                    u = queue[qhead]
                    qhead += 1
                    for v in outneighbors(dg, u)
                        if !visited[v]
                            visited[v] = true
                            qtail += 1
                            queue[qtail] = v
                        end
                    end
                end
                for k2 in 1:nms
                    psadj[k1,k2] = visited[reps[k2]]
                end
            end
        else
            error("Unknown poset_algorithm $(poset_algorithm) (expected :pairwise or :condensation)")
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

