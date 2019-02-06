include("addnoise.jl")

import ...Embeddings.AbstractEmbedding
import DelayEmbeddings.Dataset
import Simplices.Delaunay.delaunay
import StaticArrays:
    SArray, MArray

"""
    DelaunayTriangulation

Contains the vertex indices of the simplices furnishing a Delaunay triangulation.

## Indexing

If `dt` is a `DelaunayTriangulation`, then

- ``dt[1]`` returns the indices of the 1st simplex.
- ``dt[1:4]`` returns a vector of indices for the 1st through 4th simplices.
- ``dt[[5, 8]]`` returns a vector of indices for the 5th and 8th simplices.

Calling `dt[:, i]` will return the same as with single-dimension indexing, but
converts the indices into an
array where each column is a set of simplex vertex indices, i.e

``dt[:, 1]`` returns the indices of the 1st simplex.
"""
struct DelaunayTriangulation <: AbstractDelaunayTriangulation
    indices::Array{Array{Int32, 1}, 1}
end

dimension(DT::DelaunayTriangulation) = length(DT.indices[1]) - 1


"""
    DelaunayTriangulation(vertices; joggle) where {T} -> DelaunayTriangulation

Construct a Delaunay triangulation from a set of points.
"""
function DelaunayTriangulation(vertices; joggle = 0.001)
    # Slightly joggle points to avoid problems with QHull
    if joggle > 0
        addnoise!(vertices; joggle_factor = joggle)
    end
    # Get the indices of the simplices of the triangulation
    if size(vertices, 1) > size(vertices, 2)
        simplex_indices = delaunay(vertices)
    elseif size(vertices, 1) < size(vertices, 2)
        simplex_indices = delaunay(transpose(vertices))
    end

    return DelaunayTriangulation(simplex_indices)
end

function DelaunayTriangulation(vertices::Vector{AbstractVector{T}};
        joggle = 0.001) where T

    # Slightly joggle points to avoid problems with QHull
    if joggle > 0
        addnoise!(vertices; joggle_factor = joggle)
    end

    simplex_indices = delaunay(vertices)
    return DelaunayTriangulation(simplex_indices)
end

function DelaunayTriangulation(vertices::Vector{SArray{Tuple{D},T,1,D}};
        joggle = 0.001) where {D, T}

    # Slightly joggle points to avoid problems with QHull
    if joggle > 0
        addnoise!(vertices; joggle_factor = joggle)
    end

    simplex_indices = delaunay(vertices)
    return DelaunayTriangulation(simplex_indices)
end


"""
    DelaunayTriangulation(E::Embedding) -> DelaunayTriangulation

Construct a Delaunay triangulation from the points of an embedding.
"""
function DelaunayTriangulation(E::AbstractEmbedding; joggle = 0.001)
    DelaunayTriangulation(Array(transpose(E.points)), joggle = joggle)
end

"""
    DelaunayTriangulation(D::Dataset) -> DelaunayTriangulation

Construct a Delaunay triangulation from a Dataset.
"""
function DelaunayTriangulation(D::Dataset; joggle = 0.001)
    d = hcat(D.data...,)
    DelaunayTriangulation(d .+ joggle*rand(size(d)...,), joggle = joggle)
end


export DelaunayTriangulation, dimension
