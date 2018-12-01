using StateSpaceReconstruction
using Test
import DynamicalSystemsBase.Dataset

# Some example data
pts = rand(4, 25)
E = embed(pts)
D = Dataset(transpose(pts))

# Test constructors
@test DelaunayTriangulation(pts) isa DelaunayTriangulation
@test DelaunayTriangulation(transpose(pts)) isa DelaunayTriangulation
@test DelaunayTriangulation(E) isa DelaunayTriangulation
@test DelaunayTriangulation(D) isa DelaunayTriangulation

# Test indexing
DT = DelaunayTriangulation(pts)
@test DT[1] isa Vector{Int32}
@test length(DT[1]) == 5
@test DT[1:2] isa Vector{Vector{Int32}}
@test length(DT[1:2]) == 2
@test DT[:, 1] isa Array{Int32, 2}
@test size(DT[:, 1]) == (5, 1)

@test DT[:, 1:4] isa Array{Int32, 2}
@test size(DT[:, 1:4]) == (5, 4)
@test DT[:, [1, 2]] isa Array{Int32, 2}
@test DT[:, :] isa Array{Int32, 2}
