@testset "Alignment with and without zero-lagged vector matches" begin
    ts = [rand(12) for i = 1:3]
    E1 = embed(ts, [1, 2, 3, 3], [1, -1, -1, 0])
    E2 = embed(ts, [1, 2, 3], [1, -1, -1])

    @test all(E1.points[:, 1:3] .== E2.points)
    @test typeof(E1) <: AbstractEmbedding
    @test typeof(E2) <: AbstractEmbedding
end


@testset "Embedding dim > 3" begin
    ts2 = [rand(10) for i = 1:3]
    ts3 = [collect(1:10) for i = 1:5]
    # Floats
    E1 = embed(ts2, [1, 2, 3, 3], [1, -1, -1, 0])
    E2 = embed(ts2, [1, 2, 3], [1, -1, -1])

    # Integers
    E3 = embed(ts3, [1, 2, 3, 3, 1], [1, -1, -1, 0, 1])
    E4 = embed(ts3, [1, 2, 3], [1, -1, -1])

    @test all(E1.points[:, 1:3] .== E2.points)
    @test all(E3.points[:, 1:3] .== E4.points)

    @test typeof(E1) <: AbstractEmbedding
    @test typeof(E2) <: AbstractEmbedding
    @test typeof(E3) <: AbstractEmbedding
    @test typeof(E4) <: AbstractEmbedding
end

@testset "Invariantizing embeddings" begin
    pos = [1, 2, 3, 3]
    lags = [1, -1, -1, 0]
    E1 = embed([diff(rand(30)) for i = 1:4], pos, lags)
    E2 = embed([diff(rand(1:10, 30)) for i = 1:4], pos, lags)
    inv_E1 = invariantize(E1)
    inv_E2 = invariantize(E2)

    @test typeof(inv_E1) == LinearlyInvariantEmbedding{Float64}
    @test typeof(inv_E2) == LinearlyInvariantEmbedding{Int}
end


@testset "Different dispatch" begin
    u = [randn(10) for i = 1:3]
    v = [collect(1:10) for i = 1:3]
    A = randn(10, 3)
    B = hcat(v...)

    ts_inds = [1, 2, 3]
    embedding_lags = [1, 0, -2]

    # Vector of vectors
    @test typeof(embed(u)) <: AbstractEmbedding
    @test typeof(embed(v)) <: AbstractEmbedding
    @test typeof(embed(u, ts_inds, embedding_lags)) <: AbstractEmbedding
    @test typeof(embed(v, ts_inds, embedding_lags)) <: AbstractEmbedding

    # Arrays
    @test typeof(embed(A)) <: AbstractEmbedding
    @test typeof(embed(A, ts_inds, embedding_lags)) <: AbstractEmbedding

    @test typeof(embed(float.(B))) <: AbstractEmbedding
    @test typeof(embed(float.(B), ts_inds, embedding_lags)) <: AbstractEmbedding

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
