using StaticArrays
using DynamicalSystems
using RecipesBase

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
	@test typeof(E1a) <: Embeddings.AbstractEmbedding
	@test typeof(E1b) <: Embeddings.AbstractEmbedding
	@test typeof(E2a) <: Embeddings.AbstractEmbedding
	@test typeof(E2b) <: Embeddings.AbstractEmbedding
	@test typeof(E3a) <: Embeddings.AbstractEmbedding
	@test typeof(E3b) <: Embeddings.AbstractEmbedding

	#inv = invariantize(E1a, verbose = true)
	tri = delaunaytriang(E1a)
	#@test typeof(inv) <: LinearlyInvariantEmbedding
	@test typeof(tri) <: DelaunayTriangulation
end

@testset "Alignment with and without zero-lagged vector matches" begin
    ts = [rand(12) for i = 1:3]
    E1 = embed(ts, [1, 2, 3, 3], [1, -1, -1, 0])
    E2 = embed(ts, [1, 2, 3], [1, -1, -1])

    @test all(E1.points[1:3, :] .== E2.points)
    @test typeof(E1) <: Embeddings.AbstractEmbedding
    @test typeof(E2) <: Embeddings.AbstractEmbedding
end

@testset "Plotting recipes" begin
    emb = embed([diff(rand(30)) for i = 1:3], [1, 2, 3], [1, 0, -1])
    emb_invariant = invariantize(emb)

    # Test plot recipes by calling RecipesBase.apply_recipe with empty dict.
    # It should return a vector of RecipesBase.RecipeData
    d = Dict{Symbol,Any}()
    @test typeof(RecipesBase.apply_recipe(d, emb)) == Array{RecipesBase.RecipeData,1}
    @test typeof(RecipesBase.apply_recipe(d, emb_invariant)) == Array{RecipesBase.RecipeData,1}
end
