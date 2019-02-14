import StaticArrays: SVector, MVector
import DelayEmbeddings: Dataset

"""
    triangulate(pts)

Triangulate a set of points and, depending on the type of the input,
return an appropriate instance of a subtype of `AbstracTriangulationPartition`.

## Arguments
- **`pts`**. The points. May be an instance of `AbstractArray{T, 2`} (the largest dimension will be taken as the index), `Vector{Vector}`, `Vector{SVector}`, `Vector{MVector}` (each vector is a point) or a `Dataset`.
"""
function triangulate end

"""
    triangulate_m(pts)

Triangulate a set of points and, depending on the type of the input,
return an appropriate mutable instance of a subtype of `AbstractMutableTriangulationPartition`.

## Arguments
- **`pts`**. The points. May be an instance of `AbstractArray{T, 2`} (the largest dimension will be taken as the index), `Vector{Vector}`, `Vector{SVector}`, `Vector{MVector}` (each vector is a point) or a `Dataset`.
"""
function triangulate_m end


"""
    triangulate_f(pts)

Triangulate a set of points and, depending on the type of the input,
return an appropriate instance of a subtype of `AbstractTriangulationPartitionFull`,
which not only contains the points and the Delaunay triangulation, but also
the simplices and their radii, centroids, orientations and volumes.

## Arguments
- **`pts`**. The points. May be an instance of `AbstractArray{T, 2`} (the largest dimension will be taken as the index), `Vector{Vector}`, `Vector{SVector}`, `Vector{MVector}` (each vector is a point) or a `Dataset`.
"""
function triangulate_f end

"""
    triangulate_mf(pts)

Triangulate a set of points and, depending on the type of the input,
return an appropriate mutable instance of a subtype of `AbstractMutableTriangulationPartitionFull`,
which not only contains the points and the Delaunay triangulation, but also
the simplices and their radii, centroids, orientations and volumes.

## Arguments
- **`pts`**. The points. May be an instance of `AbstractArray{T, 2`} (the largest dimension will be taken as the index), `Vector{Vector}`, `Vector{SVector}`, `Vector{MVector}` (each vector is a point) or a `Dataset`.
"""
function triangulate_mf end


function triangulate(pts::Vector{Vector{T}}) where T
    DT = DelaunayTriangulation(pts)
    dim = length(pts[1])
    PointTriangulation{dim, T}(pts, DT)
end


function triangulate_m(pts::Vector{Vector{T}}) where T
    DT = DelaunayTriangulation(pts)
    dim = length(pts[1])
    MutablePointTriangulation{dim, T}(pts, DT)
end

function triangulate_f(pts::Vector{Vector{T}}) where T
    DT = DelaunayTriangulation(pts)
    n_simplices = length(DT)
    dim = length(pts[1])
    simplices = [Simplex(pts[DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    PointTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations, volumes, centroids)
end

function triangulate_mf(pts::Vector{Vector{T}}) where T
    DT = DelaunayTriangulation(pts)
    n_simplices = length(DT)
    dim = length(pts[1])
    simplices = [Simplex(pts[DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    MutablePointTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end



function triangulate(pts::AbstractArray{T, 2}) where T
    DT = DelaunayTriangulation(pts)
    dim = minimum(size(pts))

    if size(pts, 1) > size(pts, 2)
        pts = transpose(pts)
    end
    MatrixTriangulation{dim, T}(pts, DT)
end

function triangulate_m(pts::AbstractArray{T, 2}) where T
    DT = DelaunayTriangulation(pts)
    dim = minimum(size(pts))

    if size(pts, 1) > size(pts, 2)
        pts = transpose(pts)
    end

    MutableMatrixTriangulation{dim, T}(pts, DT)
end

function triangulate_f(pts::AbstractArray{T, 2}) where T
    DT = DelaunayTriangulation(pts)

    if size(pts, 1) > size(pts, 2)
        pts = transpose(pts)
    end

    n_simplices = length(DT)
    dim = minimum(size(pts))
    simplices = [Simplex(pts[:, DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    MatrixTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end

function triangulate_mf(pts::AbstractArray{T, 2}) where T
    DT = DelaunayTriangulation(pts)

    if size(pts, 1) > size(pts, 2)
        pts = transpose(pts)
    end

    n_simplices = length(DT)
    dim = minimum(size(pts))
    simplices = [Simplex(pts[:, DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    MutableMatrixTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end


function triangulate(pts::Dataset{D, T}) where {D, T}
    DT = DelaunayTriangulation(pts)
    DatasetTriangulation{D, T}(pts, DT)
end

function triangulate_m(pts::Dataset{D, T}) where {D, T}
    DT = DelaunayTriangulation(pts)
    MutableDatasetTriangulation{D, T}(pts, DT)
end

function triangulate_f(pts::Dataset{D, T}) where {D, T}
    DT = DelaunayTriangulation(pts)
    n_simplices = length(DT)
    dim = minimum(size(pts))
    simplices = [Simplex(pts[DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    MatrixTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end

function triangulate_mf(pts::Dataset{D, T}) where {D, T}
    DT = DelaunayTriangulation(pts)
    n_simplices = length(DT)
    dim = minimum(size(pts))
    simplices = [Simplex(pts[DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    MutableMatrixTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end


function triangulate(pts::Vector{SVector{D, T}}) where {D, T}
    DT = DelaunayTriangulation(pts)
    PointTriangulation{D, T}(pts, DT)
end

function triangulate_m(pts::Vector{SVector{D, T}}) where {D, T}
    DT = DelaunayTriangulation(pts)
    MutablePointTriangulation{D, T}(pts, DT)
end

function triangulate_f(pts::Vector{SVector{D, T}}) where {D, T}
    DT = DelaunayTriangulation(pts)
    n_simplices = length(DT)
    dim = length(pts[1])
    simplices = [Simplex(pts[DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    PointTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end

function triangulate_mf(pts::Vector{SVector{D, T}}) where {D, T}
    DT = DelaunayTriangulation(pts)
    n_simplices = length(DT)
    dim = length(pts[1])
    simplices = [Simplex(pts[DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    MutablePointTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end

function triangulate(pts::Vector{MVector{D, T}}) where {D, T}
    DT = DelaunayTriangulation(pts)
    PointTriangulation{D, T}(pts, DT)
end

function triangulate_m(pts::Vector{MVector{D, T}}) where {D, T}
    DT = DelaunayTriangulation(pts)
    MutablePointTriangulation{D, T}(pts, DT)
end

function triangulate_f(pts::Vector{MVector{D, T}}) where {D, T}
    DT = DelaunayTriangulation(pts)
    n_simplices = length(DT)
    dim = length(pts[1])
    simplices = [Simplex(pts[DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    PointTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end

function triangulate_mf(pts::Vector{MVector{D, T}}) where {D, T}
    DT = DelaunayTriangulation(pts)
    n_simplices = length(DT)
    dim = length(pts[1])
    simplices = [Simplex(pts[DT[i]]) for i = 1:n_simplices]
    radii = radius.(simplices)
    centroids = centroid.(simplices)
    orientations = orientation.(simplices)
    volumes = volume.(simplices)

    MutablePointTriangulationFull{dim, T}(pts, DT, simplices, radii, orientations,
        volumes, centroids)
end


export
triangulate,
triangulate_m,
triangulate_f,
triangulate_mf
