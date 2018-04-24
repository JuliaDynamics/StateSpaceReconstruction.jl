
@testset "Rectangular binning" begin
    n_pts = 10
    E1 = embed([rand(10) for i = 1:3])
    b1 = bin_equidistant(E1, n_pts)
    @test typeof(b1) == EquidistantBinning

    E2 = embed([collect(1:10) for i = 1:3])
    b2 = bin_equidistant(E2, n_pts)
    @test typeof(b2) == EquidistantBinning
end

@testset "Simplex triangulation" begin

	@testset "Triangulation" begin
        E_3D = embed([rand(30) for i = 1:3])
        E_4D = embed([rand(30) for i = 1:4])

        T_3D = triangulate(E_3D)
        @test typeof(T_3D) == Triangulation

        T_4D = triangulate(E_4D)
        @test typeof(T_4D) == Triangulation
    end

    @testset "LinearlyInvariantTriangulation" begin
        E_3D = invariantize(embed([rand(30) for i = 1:3]))
        E_4D = invariantize(embed([rand(30) for i = 1:4]))

        T_3D = triangulate(E_3D)
        @test typeof(T_3D) == LinearlyInvariantTriangulation

        T_4D = triangulate(E_4D)
        @test typeof(T_4D) == LinearlyInvariantTriangulation

    end
end
