using Test, MegaQuiverTools

@testset "Testing basic methods" begin
    K = mKroneckerquiver(4)
    
    @test MegaQuiverTools.underlying_graph(K) == [0 4;4 0]
    @test nvertices(K) == 2
    @test narrows(K) == 4
    @test isacyclic(K) == true
    @test isconnected(K) == true
    @test issource(K, 1) == true && issource(K, 2) == false
    @test issink(K, 2) == true && issink(K, 1) == false
    @test MegaQuiverTools.euler_matrix(K) == [1 -4; 0 1]
    @test Eulerform(K, [1, 1], [1, 1]) == -2
end;

@testset "Testing constructors" begin

    @test mKroneckerquiver(3).adjacency == [0 3; 0 0]
    @test threevertexquiver(3,4,5).adjacency == [0 3 4; 0 0 5; 0 0 0]
    @test loopquiver(3).adjacency == Matrix{Int64}(reshape([3],1,1))
    @test subspacequiver(2).adjacency == [0 0 1; 0 0 1;0 0 0]
end;


@testset "Testing hassemistables" begin
    Q = mKroneckerquiver(17)
    @test hassemistables(Q,[1,13], [13,-1]) == true
    @test hassemistables(Q,[1,13], [-13,1]) == false

    K = mKroneckerquiver(4)
    @test hassemistables(K,[1,5], [5,-1]) == false
    @test hassemistables(K,[1,5], [-5,1]) == false

    A2 = mKroneckerquiver(1)
    theta = [1,-1]
    d = [1,1]
    @test hassemistables(A2, d, theta) == true

    d = [2,2]
    @test hassemistables(A2, d, theta) == true

    d = [1,2]
    @test hassemistables(A2, d, theta) == false

    d = [0,0]
    @test hassemistables(A2, d, theta) == true


    K3 = mKroneckerquiver(3)
    theta = [3,-2]
    d = [2,3]
    @test hassemistables(K3, d, theta) == true

    d = [1,4]
    @test hassemistables(K3, d, theta) == false
end;

@testset "Testing allHNtypes" begin
    Q = mKroneckerquiver(3)
    d = [2,3]

    theta = [3,-2]
    @test allHNtypes(Q, d, theta) == [[[2, 3]],
    [[1, 1], [1, 2]],
    [[2, 2], [0, 1]],
    [[2, 1], [0, 2]],
    [[1, 0], [1, 3]],
    [[1, 0], [1, 2], [0, 1]],
    [[1, 0], [1, 1], [0, 2]],
    [[2, 0], [0, 3]]]

    theta = [-3,2]
    @test allHNtypes(Q, d, theta) == [[[0,3],[2,0]]]

    Q = threevertexquiver(3,4,5)
    d = [3,5,7]
    theta = [43,26,-37]
# 3vertexquiver-3-5-7-canonical.txt
    expected = "" #has to be initialised outside of the open file
    open("3vertexquiver-3-5-7-canonical.txt","r") do file
        expected = readline(file)
    end

    @test string(allHNtypes(Q,d,theta)) == expected
end;

@testset "Testing isamplystable" begin
    Q = mKroneckerquiver(3)
    d = [2,3]

    @test isamplystable(Q, d, [3,-2]) == true
    @test isamplystable(Q, d, [-3,2]) == false
end;

@testset "Testing weight handling" begin 

    U = Bundle([1,2],2)
    V = Bundle([3,3,3],3)

    @test U ⊕ V == Bundle([1, 2, 3, 3, 3], 5)
    @test U ⊗ V == Bundle([4, 4, 4, 5, 5, 5], 6)
    @test U ⊠ V == Bundle([ [1, 3], [1, 3], [1, 3],
                            [2, 3], [2, 3], [2, 3]], 6)

    @test U - V == Bundle([1, 2, -3, -3, -3], 5)

    
    @test wedge(U,0) == Bundle([0], 1)
    @test wedge(U,1) == Bundle([1, 2], 2)
    @test wedge(U,2) == Bundle([3], 1)
    @test wedge(U,3) == Bundle(Int64[], 0)
    @test wedge(U ⊠ V,2) == Bundle([[2, 6], [2, 6], [3, 6],
                                    [3, 6], [3, 6], [2, 6],
                                    [3, 6], [3, 6], [3, 6],
                                    [3, 6], [3, 6], [3, 6],
                                    [4, 6], [4, 6], [4, 6]], 15)

    @test wedge(U ⊠ V,3) == Bundle([[3, 9], [4, 9], [4, 9],
                                    [4, 9], [4, 9], [4, 9],
                                    [4, 9], [5, 9], [5, 9],
                                    [5, 9], [4, 9], [4, 9],
                                    [4, 9], [5, 9], [5, 9],
                                    [5, 9], [5, 9], [5, 9],
                                    [5, 9], [6, 9]], 20)

    @test wedge(U ⊠ V,4) == Bundle([[5, 12], [5, 12], [5, 12],
                                    [6, 12], [6, 12], [6, 12],
                                    [6, 12], [6, 12], [6, 12],
                                    [7, 12], [6, 12], [6, 12],
                                    [6, 12], [7, 12], [7, 12]], 15)

    @test wedge(U ⊠ V,5) == Bundle([[7, 15], [7, 15], [7, 15],
                                    [8, 15], [8, 15], [8, 15]], 6)

    @test wedge(U ⊠ V,6) == Bundle([[9, 18]], 1)

    W = Bundle([1, 2, 3], 3);

    @test symm(W,0) == Bundle([0], 1)
    @test symm(W,1) == Bundle([1, 2, 3], 3)
    @test symm(W,2) == Bundle([2, 3, 4, 4, 5, 6], 6)
    @test symm(W,3) == Bundle([3, 4, 5, 5, 6, 7, 6, 7, 8, 9], 10)
    @test symm(W,4) == Bundle([4, 5, 6, 6, 7, 8, 7, 8, 9, 10, 8, 9, 10, 11, 12], 15)
end;