export isolate_cellsubsets

"""
    isolate_cellsubsets(ec::EuclideanComplex, mvf::CellSubsets, subsets::CellSubsets)

Isolate the dynamics associated with a selected collection of Morse sets.

Given a Euclidean complex `ec`, a multivector field `mvf`, and a collection of
Morse sets `subsets`, this function computes the smallest isolated invariant set
containing the selected Morse sets and their connecting orbits, fattens it by a
delta neighborhood to ensure local closedness, and returns the restricted
complex and induced multivector field.

The input `subsets` should contain Morse sets as determined by `morse_sets`.
The minimum pairwise distance delta between the selected Morse sets is used to
fatten the isolated invariant set by including all 2-cells incident to a
vertex within delta of it.

Returns `(ec_res, mvf_res)`, the restricted Euclidean complex and multivector
field capturing the dynamics between the selected Morse sets.

A warning is issued if the fattened set is not a single coherent multivector,
which may indicate that the set is not locally closed.

# Example
```julia
ec_res, mvf_res = isolate_cellsubsets(ecr, mvfr, morsedecompr[[1, 3, 5]])
morsedecomp_res, mdposet_res = morse_sets(ec_res, mvf_res, poset=true)
```
"""
function isolate_cellsubsets(ec::EuclideanComplex,
                             mvf::CellSubsets,
                             subsets::CellSubsets)

    n = length(subsets)

    if n < 2
        error("At least two Morse sets are required for isolation.")
    end

    ##########################################################################
    # Step 1: Compute the smallest isolated invariant set (Morse interval)   #
    # containing all selected Morse sets and their connecting orbits         #
    ##########################################################################

    u_set = morse_interval(ec, mvf, subsets)

    ##########################################################################
    # Step 2: Compute minimum pairwise distance delta between selected sets  #
    ##########################################################################

    # 0-skeleton (vertex indices) for each selected Morse set

    skel = [lefschetz_skeleton(ec, subsets[i], 0) for i in 1:n]

    # Compute minimum pairwise distance between selected Morse sets

    deez = Float64[]
    for i = 1:n-1
        for j = 1:n-i
            d = Float64[]
            for v in skel[i+j]
                dpoint = ec.coords[v][1]
                dv = cellsubset_distance(ec, subsets[i], dpoint)
                if dv > 0
                    push!(d, dv)
                end
            end
            # Skip pairs where every vertex of one set coincides with the
            # other (e.g. two directly touching Morse sets) — they impose
            # no lower bound on delta, but other pairs still can.
            if !isempty(d)
                push!(deez, minimum(d))
            end
        end
    end

    if isempty(deez)
        error("All selected Morse sets are mutually touching (zero distance) — " *
              "cannot determine a fattening delta.")
    end

    delta = minimum(deez)
    println("Delta: $delta")

    ##########################################################################
    # Step 3: Fatten u_set by adding 2-cells within delta                   #
    ##########################################################################

    # 0-skeleton of all 2-cells in ec

    faces = findall(isequal(2), ec.dimensions)
    f_skel = lefschetz_skeleton(ec, faces, 0)

    # Add all triangles incident to a vertex within delta of u_set

    B = Vector{Int64}[]
    for v in f_skel
        dpoint = ec.coords[v][1]
        if cellsubset_distance(ec, u_set, dpoint) <= delta
            cell = lefschetz_openhull(ec, [v])
            cell = lefschetz_skeleton(ec, cell, 2)
            push!(B, cell)
        end
    end

    u_fat = isempty(B) ? u_set : union(reduce(union, B), u_set)

    ##########################################################################
    # Step 4: Create, extract, and restrict                                  #
    ##########################################################################

    A = deepcopy(mvf)
    S = push!(A, u_fat)

    mvs = create_mvf_hull(ec, S)
    mvu = extract_multivectors(ec, mvs, u_fat)

    if length(mvu) != 1
        @warn "u_fat is not a single coherent multivector — it may not be locally closed. Proceeding anyway."
    end

    cells = mvu[1]

    ec_res, mvf_res = restrict_dynamics(ec, mvf, cells)

    return ec_res, mvf_res

end
