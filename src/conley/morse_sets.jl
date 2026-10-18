export morse_sets

"""
    morse_sets(lc::AbstractComplex, mvf::CellSubsets; poset::Bool=false,
               poset_algorithm::Symbol=:scc_firsthit)

Find the nontrivial Morse sets of a multivector field on a Lefschetz complex.

The input argument `lc` contains the Lefschetz complex, and `mvf` describes
the multivector field. The function returns the nontrivial Morse sets as
a `Vector{Vector{Int}}`. If the optional argument `poset=true` is added,
then the function returns both the Morse sets and the adjacency matrix of
the Hasse diagram of the underlying poset.

`poset_algorithm` selects how the Hasse diagram is computed when
`poset=true`:
* `:scc_firsthit` (default) computes the Hasse diagram directly in a
  single topological-order pass over the Tarjan condensation of the
  cell-level digraph (one node per strongly connected component), using
  a first-hit dynamic-programming recurrence -- no `nms^2` reachability
  matrix is ever materialized, and no separate `transitivereduction`
  post-processing step is needed. This is the fastest and current
  upstream-default algorithm (from ConleyDynamics 0.7.25).
* `:pairwise` issues one `has_path` query per ordered pair of Morse sets
  (`nms^2` total), each a fresh traversal of the full cell-level
  digraph, followed by a separate `transitivereduction` call. Simple,
  well-tested, but can dominate runtime once `nms` grows into the
  thousands -- kept for benchmarking/comparison against the newer
  algorithms.
* `:condensation` computes the same `nms^2` adjacency matrix as
  `:pairwise`, but via one single-source reachability traversal per
  Morse set (`nms` total, parallelized across threads) instead of
  `nms^2` -- the reachable set from each representative is computed
  once and reused for every membership check against the other
  representatives, then fed through the same `transitivereduction` step
  as `:pairwise`. Produces identical results to `:pairwise`, cheaper
  for large `nms`, but still pays for the reachability matrix and the
  transitive-reduction step that `:scc_firsthit` avoids entirely. Kept
  for benchmarking/comparison.
"""
function morse_sets(lc::AbstractComplex, mvf::CellSubsets; poset::Bool=false,
                    poset_algorithm::Symbol=:scc_firsthit)
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
    # can be determined from the multivector cells alone. The
    # algorithm "auto" uses "rank" for small multivectors and
    # "morse" for large ones.

    Threads.@threads for k = 1:length(scc)
        if length(unique(cell2mv[scc[k]])) > 1
            isnontrivial[k] = 1
        else
            isnontrivial[k] = sum(conley_index(lc, scc[k], algorithm="auto",
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

        if poset_algorithm == :pairwise
            psadj = fill(false,nms,nms)
            for k1 in 1:nms
                vk1 = reps[k1]
                for k2 in 1:nms
                    vk2 = reps[k2]
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
        elseif poset_algorithm == :condensation
            # One single-source reachability traversal per Morse set
            # (parallelized) instead of nms^2 individual has_path calls --
            # the reachable set from each representative is computed once
            # and reused for every membership check, instead of re-deriving
            # it from scratch for every (k1,k2) pair. Manual BFS with an
            # index-based queue (not popfirst!) to avoid O(n)-per-dequeue
            # overhead on a plain Vector.
            nv_dg = nv(dg)
            psadj = fill(false,nms,nms)
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
            psdg = SimpleDiGraph(psadj)
            pshasse = transitivereduction(psdg)

            hasse = fill(false,nms,nms)
            for k1 in 1:nms
                for k2 in 1:nms
                    hasse[k2,k1] = has_edge(pshasse, k1, k2)
                end
            end
        elseif poset_algorithm == :scc_firsthit
            # Advisor's fix (upstream 0.7.25): compute the Hasse diagram
            # directly in a single topological-order pass over the
            # condensation (one node per strongly connected component),
            # instead of ever materializing an nms x nms reachability
            # matrix and then transitively reducing it. For each
            # component, "firsthit" collects the Morse sets reachable
            # along paths that do not pass through any other Morse set
            # first; a first-hit Morse set survives into the Hasse
            # diagram unless it is also reachable via another first-hit
            # Morse set (the `covered` check), computed with a `stamp`
            # array to dedupe targets in O(1) instead of a Set per node.

            nscc = length(scc)

            # Map the cells to the strongly connected components,
            # and the components to the Morse sets (zero if none)

            cell2scc = fill(Int(0), nv(dg))
            for k = 1:nscc
                for v in scc[k]
                    cell2scc[v] = k
                end
            end

            scc2ms = fill(Int(0), nscc)
            scc2ms[ntindices] = 1:nms

            # Process the components in reverse topological order
            # of the condensation, i.e., successors come first

            sccorder = reverse(topological_sort(condensation(dg, scc)))

            # For every component, find the Morse sets which can be reached
            # along paths that do not pass through any other Morse set

            firsthit = [Int[] for _ in 1:nscc]
            stamp    = fill(Int(0), nms)

            for k in sccorder
                fh = firsthit[k]
                for v in scc[k]
                    for w in outneighbors(dg, v)
                        kw = cell2scc[w]
                        kw == k && continue
                        targets = scc2ms[kw] > 0 ? (scc2ms[kw]:scc2ms[kw]) : firsthit[kw]
                        for m in targets
                            if stamp[m] != k
                                stamp[m] = k
                                push!(fh, m)
                            end
                        end
                    end
                end
            end

            # The transitive closure of the first-hit relation is the order
            # between the Morse sets. A first-hit Morse set is in the Hasse
            # diagram unless it can be reached via another first-hit Morse set

            reach   = [falses(nms) for _ in 1:nms]
            covered = falses(nms)
            hasse   = fill(false,nms,nms)

            for k in sccorder
                m1 = scc2ms[k]
                m1 == 0 && continue
                fill!(covered, false)
                for m2 in firsthit[k]
                    covered .|= reach[m2]
                end
                for m2 in firsthit[k]
                    reach[m1][m2] = true
                    if !covered[m2]
                        hasse[m2,m1] = true
                    end
                end
                reach[m1] .|= covered
            end
        else
            error("Unknown poset_algorithm $(poset_algorithm) (expected :scc_firsthit, :pairwise, or :condensation)")
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
