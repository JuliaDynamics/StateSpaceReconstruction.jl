
@testset "Rectangular binning" begin
    n_pts = 30
    n_bins = 4
    E1 = embed([rand(n_pts) for i = 1:3])
    b1 = bin_equidistant(E1, n_bins) # bin sizes
    b2 = bin_equidistant(E1, 0.1) # bin sizes as fraction of ranges along each coord. axis
    b3 = bin_equidistant(E1, [0.1, 0.1, 0.2])

    @test typeof(b1) == EquidistantBinning
    @test typeof(b2) == EquidistantBinning
    @test typeof(b3) == EquidistantBinning

    E2 = embed([collect(1:10) for i = 1:3])
    b2 = bin_equidistant(E2, n_pts)
    @test typeof(b2) == EquidistantBinning
end

@testset "Simplex triangulation" begin
    n_pts = 30
	@testset "Triangulation" begin
        E_3D = embed([rand(n_pts) for i = 1:3])
        E_4D = embed([rand(n_pts) for i = 1:4])

        T_3D = triangulate(E_3D)
        @test typeof(T_3D) == GenericTriangulation

        T_4D = triangulate(E_4D)
        @test typeof(T_4D) == GenericTriangulation
    end

    @testset "LinearlyInvariantTriangulation" begin
        E_3D = invariantize(embed([rand(n_pts) for i = 1:3]))
        E_4D = invariantize(embed([rand(n_pts) for i = 1:4]))

        T_3D = triangulate(E_3D)
        @test typeof(T_3D) == LinearlyInvariantTriangulation

        T_4D = triangulate(E_4D)
        @test typeof(T_4D) == LinearlyInvariantTriangulation

    end
end
