using ConleyDynamics
using Test

@testset "Tutorial" begin
    #
    # This test set covers the tutorial examples
    #
    labels = ["A","B","C","D","E","F"]
    simplices = [["A","B"],["A","C"],["B","C"],["B","D"],["D","E","F"]]
    sc = create_simplicial_complex(labels,simplices)

    @test length(sc.labels) == 14
    @test sum(sc.dimensions) == 9
    @test sparse_size(sc.boundary,1) == sparse_size(sc.boundary,2)
    @test sparse_size(sc.boundary,1) == 14
    @test sc.ncells == 14
    @test sc.dim == 2

    @test homology(sc) == [1,1,0]
    @test relative_homology(sc, [1,6]) == [0,2,0]
    @test relative_homology(sc, ["DE","DF","EF"]) == [0,1,1]
    @test relative_homology(sc, []) == [1,1,0]

    filtration = [1,1,1,2,2,2,1,1,1,3,2,2,2,4]
    phsingles, phpairs = persistent_homology(sc, filtration)

    @test phsingles == [[1],[1],[]]
    @test phpairs == [[(2, 3)], [(2, 4)], []]

    formanvf = [["A","AC"],["B","AB"],["C","BC"],["D","BD"],["E","DE"]]
    formaninfo = mvf_information(sc,formanvf)

    @test formaninfo["N mv"] == 9
    @test formaninfo["N regular"] == 5
    @test formaninfo["Lengths regular"] == [[2, 5]]
    @test formaninfo["Lengths critical"] == [[1, 4]]
    @test formaninfo["N critical"] == 4
    @test conley_index(sc, ["F"]) == [1,0,0]
    @test conley_index(sc, ["DF"]) == [0,1,0]
    @test conley_index(sc, ["DEF"]) == [0,0,1]
    @test conley_index(sc, ["AB", "AC", "BC", "A", "B", "C"]) == [1,1,0]

    cm = connection_matrix(sc, formanvf)

    @test cm.morse == [["A","B","C","AB","AC","BC"],["F"],["DF"],["EF"],["DEF"]]
    @test cm.conley == [[1,1,0],[1,0,0],[0,1,0],[0,1,0],[0,0,1]]
    @test sum(full_from_sparse(cm.matrix)) == 6

    labels = ["A","B","C","D"]
    simplices = [["A","B","C"],["B","C","D"]]
    sclogo = create_simplicial_complex(labels,simplices)
    mvflogo = [["A","AB"],["C","AC"],["B","BC","BD","BCD"]]
    
    cl1, mo1 = lefschetz_clomo_pair(sclogo, ["ABC"])
    cl2, mo2 = lefschetz_clomo_pair(sclogo, mvflogo[3])
    cmlogo = connection_matrix(sclogo, mvflogo)

    @test relative_homology(sclogo, cl1, mo1) == [0,0,1]
    @test conley_index(sclogo, ["ABC"]) == [0,0,1]
    @test relative_homology(sclogo, cl2, mo2) == [0,0,0]
    @test cmlogo.conley == [[1,0,0],[0,1,0],[0,0,1]]
    @test full_from_sparse(cmlogo.matrix) == [0 0 0;0 0 1;0 0 0]
end





