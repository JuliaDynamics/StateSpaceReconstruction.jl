using DynamicalSystems.Dataset

@testset "Alignment with and without zero-lagged vector matches" begin
    ts = [rand(12) for i = 1:3]
    E1 = embed(ts, [1, 2, 3, 3], [1, -1, -1, 0])
    E2 = embed(ts, [1, 2, 3], [1, -1, -1])

    @test all(E1.points[:, 1:3] .== E2.points)
    @test typeof(E1) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(E2) <: StateSpaceReconstruction.AbstractEmbedding
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

    @test typeof(E1) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(E2) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(E3) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(E4) <: StateSpaceReconstruction.AbstractEmbedding
end


@testset "Multiple dispatch" begin
    u = [randn(10) for i = 1:3]
    v = [collect(1:10) for i = 1:3]
    A = randn(10, 3)
    B = hcat(v...)
    D = Dataset(A)
    ts_inds = [1, 2, 3]
    embedding_lags = [1, 0, -2]

    # Vector of vectors
    @test typeof(embed(u)) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(embed(v)) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(embed(u, ts_inds, embedding_lags)) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(embed(v, ts_inds, embedding_lags)) <: StateSpaceReconstruction.AbstractEmbedding

    # Arrays
    @test typeof(embed(A)) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(embed(A, ts_inds, embedding_lags)) <: StateSpaceReconstruction.AbstractEmbedding

    @test typeof(embed(float.(B))) <: StateSpaceReconstruction.AbstractEmbedding
    @test typeof(embed(float.(B), ts_inds, embedding_lags)) <: StateSpaceReconstruction.AbstractEmbedding

    # Datasets
    @test typeof(embed(D)) <: Embeddings.AbstractEmbedding
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
