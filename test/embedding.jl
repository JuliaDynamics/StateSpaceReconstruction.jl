ts = [SingleTimeSeries(rand(10)) for i = 1:3]
ts2 = [rand(10) for i = 1:3]
ts3 = [collect(1:10) for i = 1:5]
@testset "Alignment with and without zero-lagged vector matches" begin
    E1 = embed(ts, [1, 2, 3, 3], [1, -1, -1, 0])
    E2 = embed(ts, [1, 2, 3], [1, -1, -1])

    @test all(E1.points[:, 1:3] .== E2.points)
end


@testset "Embedding dim > 3" begin

    # Floats
    E1 = embed(ts2, [1, 2, 3, 3], [1, -1, -1, 0])
    E2 = embed(ts2, [1, 2, 3], [1, -1, -1])

    @test all(E1.points[:, 1:3] .== E2.points)

    # Ints
    E3 = embed(ts3, [1, 2, 3, 3, 1], [1, -1, -1, 0, 1])
    E4 = embed(ts3, [1, 2, 3], [1, -1, -1])

    @test all(E3.points[:, 1:3] .== E4.points)
end
