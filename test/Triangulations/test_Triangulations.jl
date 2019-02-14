using StateSpaceReconstruction
import DelayEmbeddings.Dataset

import StateSpaceReconstruction:
    MutablePointTriangulation, PointTriangulation,
    MutableDatasetTriangulation, DatasetTriangulation,
    MutableMatrixTriangulation, MatrixTriangulation,
    DelaunayTriangulation,
    triangulate, triangulate_m

import StaticArrays: SVector, MVector

using DelayEmbeddings

pts = [rand(3) for i = 1:30];

@testset "Triangulation structs" begin

    @test triangulate(Dataset(pts)) isa DatasetTriangulation{3,Float64}
    @test triangulate([SVector{3, Float64}(x) for x in pts]) isa PointTriangulation{3, Float64}
    @test triangulate([MVector{3, Float64}(x) for x in pts]) isa PointTriangulation{3, Float64}
    @test triangulate(hcat(pts...,)) isa MatrixTriangulation{3,Float64}
    @test triangulate(transpose(hcat(pts...,))) isa MatrixTriangulation{3,Float64}

    @test triangulate_m(Dataset(pts)) isa MutableDatasetTriangulation{3,Float64}
    @test triangulate_m([SVector{3, Float64}(x) for x in pts]) isa MutablePointTriangulation{3, Float64}
    @test triangulate_m([MVector{3, Float64}(x) for x in pts]) isa MutablePointTriangulation{3, Float64}
    @test triangulate_m(hcat(pts...,)) isa MutableMatrixTriangulation{3,Float64}
    @test triangulate_m(transpose(hcat(pts...,))) isa MutableMatrixTriangulation{3,Float64}
end
