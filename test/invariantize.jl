
@testset "Invariantizing embeddings" begin
    pos = [1, 2, 3, 3]
    lags = [1, -1, -1, 0]
    E1 = embed([diff(rand(30)) for i = 1:4], pos, lags)
    E2 = embed([diff(rand(1:10, 30)) for i = 1:4], pos, lags)
    inv_E1 = invariantize(E1)
    @test typeof(inv_E1) == LinearlyInvariantEmbedding{4, Float64}
end


@testset "Invariantizing points" begin
    pos = [1, 2, 3, 3]
    lags = [1, -1, -1, 0]
    pts = rand(3, 30)
    inv_pts = invariantize(pts)
    @test forwardlinearmap_invariant(inv_pts) == true
end

# Use Simplices.SimplexSplitting until the splitting routines are re-implemented
# for the column-wise approach
using Simplices.SimplexSplitting

"""
    canonical_simplex(dim; split_factor = 1)

Returns a tuple containing the vertices
"""
function canonical_simplex(dim; split_factor = 1)
    canonical_simplex_vertices = zeros(E + 1, E)
    canonical_simplex_vertices[2:(E+1), :] = Matrix(1.0I, E, E)
    simplex_indices = zeros(Int, 1, E + 1)
    simplex_indices[1, :] = round.(Int, collect(1:E+1))

    refined = refine_triangulation(canonical_simplex_vertices, simplex_indices, [1], k)
    triang_vertices, triang_simplex_indices = refined[1], refined[2]
end
