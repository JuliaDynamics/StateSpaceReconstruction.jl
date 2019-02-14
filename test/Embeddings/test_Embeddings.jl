using StaticArrays
using Distances
using NearestNeighbors

function generate_embeddings(n, dim)
	A = rand(n, dim)
	S = SMatrix{n, dim}(A);
	v = [rand(n) for i = 1:dim]
	sv = [S[:, i] for i = 1:dim]
	positions = [rand(1:dim) for i = 1:dim]
	lags = [i for i = 1:dim]

	# Regular arrays
	E1a = cembed(A, positions, lags)
	E1b = cembed(A)

	# Static arrays
	E2a = cembed(S, positions, lags)
	E2b = cembed(S)


	E1a, E1b, E2a, E2b
end

@testset "Embedding in $D-D" for D = 2:4

	E1a, E1b, E2a, E2b = generate_embeddings(50, D)
	@test typeof(E1a) <: Embeddings.AbstractEmbedding
	@test typeof(E1b) <: Embeddings.AbstractEmbedding
	@test typeof(E2a) <: Embeddings.AbstractEmbedding
	@test typeof(E2b) <: Embeddings.AbstractEmbedding

	#inv = invariantize(E1a, verbose = true)
	#tri = DelaunayTriangulation(E1a)
	#@test typeof(inv) <: LinearlyInvariantEmbedding
	#@test typeof(tri) <: DelaunayTriangulation
end

@testset "Alignment with and without zero-lagged vector matches" begin
    ts = [rand(12) for i = 1:3]
    E1 = cembed(ts, [1, 2, 3, 3], [1, -1, -1, 0])
    E2 = cembed(ts, [1, 2, 3], [1, -1, -1])

    @test all(E1.points[1:3, :] .== E2.points)
    @test typeof(E1) <: Embeddings.AbstractEmbedding
    @test typeof(E2) <: Embeddings.AbstractEmbedding
end


@testset "Interface with NearestNeighbors.jl" begin
	E = StateSpaceReconstruction.cembed([rand(100)], [1, 1, 1, 1], [0, -3, -4, 5])
	brutetree = BruteTree(E, Euclidean())
	balltree = BallTree(E, Euclidean())
	kdtree = KDTree(E, Minkowski(2.5))
	@test typeof(brutetree) <: NearestNeighbors.BruteTree
	@test typeof(balltree) <: NearestNeighbors.BallTree
	@test typeof(kdtree) <: NearestNeighbors.KDTree
end
