@testset "Topological operations" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # Closure of a single edge includes its two vertices
    cl_AB = lefschetz_closure(sc, ["AB"])
    @test sort(cl_AB) == sort(["A","B","AB"])

    # Closure of a triangle includes all faces
    cl_DEF = lefschetz_closure(sc, ["DEF"])
    @test sort(cl_DEF) == sort(["D","E","F","DE","DF","EF","DEF"])

    # Closure of a vertex is just the vertex
    @test lefschetz_closure(sc, ["A"]) == ["A"]

    # Skeleton: dimension 0 returns only vertices (as integer indices)
    sk0 = lefschetz_skeleton(sc, 0)
    @test length(sk0) == 6
    @test all(sc.dimensions[v] == 0 for v in sk0)

    # Skeleton: dimension 1 returns only edges, not vertices (as integer indices)
    sk1 = lefschetz_skeleton(sc, 1)
    @test all(sc.dimensions[v] == 1 for v in sk1)
    @test length(sk1) == 7
    @test !(sc.indices["A"] in sk1)

    # Neighbors of vertex F: cells in closure or openhull of {F} but not F itself
    nb_F = lefschetz_neighbors(sc, ["F"])
    @test "DF" in nb_F
    @test "EF" in nb_F
    @test "DEF" in nb_F
    @test !("F" in nb_F)
end

@testset "Open and closed predicates" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # Closure is closed
    @test lefschetz_is_closed(sc, lefschetz_closure(sc, ["AB"]))

    # Single edge without vertices is not closed
    @test !lefschetz_is_closed(sc, ["AB"])

    # Vertices alone are closed
    @test lefschetz_is_closed(sc, ["A","B","C"])

    # Full complex is closed
    @test lefschetz_is_closed(sc, sc.labels)

    # Empty set is closed
    @test lefschetz_is_closed(sc, String[])

    # A closed set is always locally closed
    @test lefschetz_is_locally_closed(sc, lefschetz_closure(sc, ["DEF"]))

    # Single edge (without vertices) is locally closed (mouth = {A,B} is closed)
    @test lefschetz_is_locally_closed(sc, ["AB"])
end

@testset "Subcomplex extraction" begin
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels, simplices)

    # Closed subcomplex of triangle DEF: 7 cells
    cl_DEF = lefschetz_closure(sc, ["DEF"])
    sub = lefschetz_subcomplex(sc, cl_DEF)
    @test sub.ncells == 7
    # Filled triangle is contractible: [1, 0, 0]
    @test homology(sub) == [1, 0, 0]

    # Extract the boundary circle of the triangle: 3 vertices + 3 edges
    circle = setdiff(cl_DEF, ["DEF"])
    sub2 = lefschetz_subcomplex(sc, circle)
    @test sub2.ncells == 6
    # Circle has homology [1, 1]
    @test homology(sub2) == [1, 1]

    # The subcomplex boundary is the minor of the full boundary matrix,
    # also over the rationals, and validation can be skipped
    sc0 = create_simplicial_complex(labels, simplices, p=0)
    for set in (circle, cl_DEF, ["DE","DEF","EF"], ["AB","AC","BC","A","B","C"])
        I = sort(convert_cells(sc0, set))
        sub0 = lefschetz_subcomplex(sc0, set)
        @test full_from_sparse(sub0.boundary) ==
              full_from_sparse(sparse_minor(sc0.boundary, I, I))
        @test sub0.labels == sc0.labels[I]
        sub0n = lefschetz_subcomplex(sc0, set, validate=false)
        @test full_from_sparse(sub0n.boundary) == full_from_sparse(sub0.boundary)
    end
    @test_throws ErrorException lefschetz_subcomplex(sc, ["D","DEF"])
end

@testset "Surfaces" begin
    # Torus: H = [1, 2, 1] over any field
    @test homology(simplicial_torus(0)) == [1, 2, 1]
    @test homology(simplicial_torus(2)) == [1, 2, 1]

    # Klein bottle: H = [1, 1, 0] over rationals, [1, 2, 1] over GF(2)
    @test homology(simplicial_klein_bottle(0)) == [1, 1, 0]
    @test homology(simplicial_klein_bottle(2)) == [1, 2, 1]

    # Projective plane: H = [1, 0, 0] over rationals, [1, 1, 1] over GF(2)
    @test homology(simplicial_projective_plane(0)) == [1, 0, 0]
    @test homology(simplicial_projective_plane(2)) == [1, 1, 1]
end

@testset "Interior, filtrations, and reductions" begin
    for lc in (create_cubical_rectangle(4, 3)[1], create_simplicial_rectangle(3, 3)[1],
               example_clorenz()[1])
        N = lc.ncells
        cl = [lefschetz_closure(lc, [d]) for d in 1:N]
        # Deterministic pseudo-random subsets and filtration values
        pseudo(k, m) = mod(7919 * k + 104729 * m, 101)

        # Interior: cells whose star is contained in the subset
        for m in 1:3
            sub  = [c for c in 1:N if pseudo(c, m) < 30 * m]
            ref  = [c for c in sub if all(d -> !(c in cl[d]) || (d in sub), 1:N)]
            @test sort(lefschetz_interior(lc, sub)) == ref
            @test sort(lefschetz_interior(lc, lc.labels[sub])) == sort(lc.labels[ref])
            @test sort(lefschetz_topboundary(lc, sub)) ==
                  setdiff(lefschetz_closure(lc, sub), ref)
        end

        # Filtration: every cell enters at the smallest positive value
        # of a cell whose closure contains it (0 if there is none)
        for (lo, K) in ((0, 3), (1, 4), (1, N ÷ 2))
            fv = [lo + mod(pseudo(c, K), K - lo + 1) for c in 1:N]
            fv[1:K-lo+1] = lo:K
            newf = [minimum((fv[d] for d in 1:N if (c in cl[d]) && (fv[d] > 0)); init=typemax(Int))
                    for c in 1:N]
            newf[newf .== typemax(Int)] .= 0
            lcsub, fvsub = lefschetz_filtration(lc, fv)
            @test lcsub.labels == lc.labels[newf .> 0]
            @test fvsub == newf[newf .> 0]
        end

        # Persistent homology with an injective filter
        filt = create_random_filter(lc)
        ph = persistent_homology(lc, filt, algorithm="matrix")
        @test persistent_homology(lc, filt, algorithm="morse") == ph
        @test persistent_homology(lc, filt, algorithm="pmorse") == ph

        # Reductions along disjoint vertex-edge pairs preserve homology
        used  = fill(false, N)
        pairs = Vector{Vector{Int}}()
        for e in 1:N
            if lc.dimensions[e] == 1
                vs = lc.boundary.columns[e]
                if !any(used[vs])
                    push!(pairs, isodd(e) ? [vs[1], e] : [e, vs[1]])
                    used[vs] .= true
                end
            end
        end
        removed = reduce(vcat, pairs)
        lcred = lefschetz_reduction(lc, pairs)
        @test lcred.labels == lc.labels[setdiff(1:N, removed)]
        @test homology(lcred) == homology(lc)
        lcredL = lefschetz_reduction(lc, [lc.labels[p] for p in pairs])
        @test full_from_sparse(lcredL.boundary) == full_from_sparse(lcred.boundary)
    end
end
