using Reexport
@reexport module Delaunay

import Simplices.Delaunay.delaunay
import ...Embeddings
using Statistics
using Distributions
import DynamicalSystemsBase.Dataset

####################################
# Abstract type
####################################


abstract type AbstractDelaunayTriangulation end

ADT = AbstractDelaunayTriangulation
dimension(DT::ADT) = length(DT.indices[1]) - 1
nsimplices(DT::ADT) = length(DT.indices)
Base.size(DT::ADT) = (dimension(DT), nsimplices(DT))
Base.length(DT::ADT) = length(DT.indices)

# Indexing
Base.getindex(DT::ADT, i) = DT.indices[i]
Base.getindex(DT::ADT, i::Colon, j::Colon) = hcat(DT.indices...,)
Base.getindex(DT::ADT, i::Colon, j) = hcat(DT[j]...,)
Base.getindex(DT::ADT, i::Colon, j::Int) = hcat(DT[j])

Base.firstindex(DT::ADT) = 1
Base.lastindex(DT::ADT) = length(DT)

function summarise(DT::AbstractDelaunayTriangulation)
    _type = typeof(DT)
    n_simplices = nsimplices(DT)
    D = dimension(DT)
    summary = "$D-dimensional $_type with $n_simplices simplices\n"
end

Base.show(io::IO, DT::AbstractDelaunayTriangulation) = println(io, summarise(DT))

@inline Base.length(DT::AbstractDelaunayTriangulation) = nsimplices(DT)



####################################
# Composite types
####################################

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

"""
    addnoise!(pts, joggle_factor)

Adding uniformly distributed noise to each observation equal to `joggle_factor`
times the maximum of the standard deviations for each variables.
"""
function addnoise!(pts; joggle_factor = 0.001)

    D = minimum(size(pts))
    npts = maximum(size(pts))

    if size(pts, 1) > size(pts, 2)
        dim = 1
    else
        dim = 2
    end

    # Scale standard deviation along each axis by joggle factor
    σ = joggle_factor .* std(pts, dims = dim)

    for i in 1:D
        r = [rand(Uniform(-σ[i], σ[i])) for pt in 1:npts]
        if dim == 1
            pts[:, i] .+= r
        elseif dim == 2
            pts[i, :] .+= r
        end
    end
end


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


"""
    DelaunayTriangulation(E::Embedding) -> DelaunayTriangulation

Construct a Delaunay triangulation from the points of an embedding.
"""
function DelaunayTriangulation(E::Embeddings.AbstractEmbedding; joggle = 0.001)
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


export
AbstractDelaunayTriangulation,
DelaunayTriangulation,
nsimplices,
dimension

end
