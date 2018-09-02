using StaticArrays, DynamicalSystems

function generate_embeddings(n, dim)
	A = rand(n, dim)
	S = SMatrix{n, dim}(A);
	D = Dataset(Array(S));
	v = [rand(n) for i = 1:dim]
	sv = [S[:, i] for i = 1:dim]
	positions = [rand(1:dim) for i = 1:dim]
	lags = [i for i = 1:dim]

	# Regular arrays
	E1a = embed(A, positions, lags)
	E1b = embed(A)

	# Static arrays
	E2a = embed(S, positions, lags)
	E2b = embed(S)

	# DynamicalSystemsBase.Dataset
	E3a = embed(D, positions, lags)
	E3b = embed(D)

	E1a, E1b, E2a, E2b, E3a, E3b
end

@testset "Embedding in $D-D" for D = 2:4

	E1a, E1b, E2a, E2b, E3a, E3b = generate_embeddings(50, D)
	@test typeof(E1a) <: AbstractEmbedding
	@test typeof(E1b) <: AbstractEmbedding
	@test typeof(E2a) <: AbstractEmbedding
	@test typeof(E2b) <: AbstractEmbedding
	@test typeof(E3a) <: AbstractEmbedding
	@test typeof(E3b) <: AbstractEmbedding

	inv = invariantize(E1a, verbose = true)
	tri = delaunaytriang(E1a)
	@test typeof(inv) <: LinearlyInvariantEmbedding
	@test typeof(tri) <: DelaunayTriangulation


end
