@testset "Assign bin labels" begin
    D = 3
    E = embed([rand(30) for i = 1:D])
    npts = size(E.points, 1)

    @testset "ϵ is an Int" begin
        labels = assign_bin_labels(E::AbstractEmbedding, 10)
        @test size(labels, 1) == D
        @test size(labels, 2) == npts
    end

    @testset "ϵ is a Float64" begin
        labels = assign_bin_labels(E::AbstractEmbedding, 0.2)
        @test size(labels, 1) == D
        @test size(labels, 2) == npts
    end

    @testset "ϵ is a Vector{Float64}" begin
        labels = assign_bin_labels(E::AbstractEmbedding, [0.2, 0.3, 0.1])
        @test size(labels, 1) == D
        @test size(labels, 2) == npts
    end

    @testset "ϵ is a Vector{Int}" begin
        labels = assign_bin_labels(E::AbstractEmbedding, [3, 4, 2])
        @test size(labels, 1) == D
        @test size(labels, 2) == npts
    end
end

@testset "Rectangular binning" begin
    n_pts = 30
    n_bins = 4
    E1 = embed([rand(n_pts) for i = 1:3])
    b1 = bin_rectangular(E1, n_bins) # bin sizes
    b2 = bin_rectangular(E1, 0.1) # bin sizes as fraction of ranges along each coord. axis
    b3 = bin_rectangular(E1, [0.1, 0.1, 0.2])

    @test typeof(b1) <: RectangularBinning
    @test typeof(b2) <: RectangularBinning
    @test typeof(b3) <: RectangularBinning

    E2 = embed([collect(1:10) for i = 1:3])
    b2 = bin_rectangular(E2, n_pts)
    @test typeof(b2) <: RectangularBinning
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
