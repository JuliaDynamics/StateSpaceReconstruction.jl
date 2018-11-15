@reexport module Delaunay

import Simplices.Delaunay: delaunay
using ..Embeddings
using StaticArrays: SVector
using DynamicalSystemsBase: Dataset

"""
    DelaunayTriangulation{D, T}

A Delaunay triangulation in dimension D. If `d`
is an instance of `DelaunayTriangulation`, then
`indices[i]` gives the D + 1 indices of the vertices
corresponding to the i-th simplex. The indices are
expressed in terms of the points it was produced
from.
"""
struct DelaunayTriangulation{D, T}
    indices::Dataset{D, T}
end

@inline Base.length(d::DelaunayTriangulation{D,T}) where {D,T} = length(d.indices)
@inline Base.size(d::DelaunayTriangulation{D,T}) where {D,T} = (length(d), D)
@inline Base.size(d::DelaunayTriangulation, i) = size(d.indices)[i]

@inline Base.eachindex(d::DelaunayTriangulation) = Base.OneTo(length(d.indices))

tps = Union{SVector{D, T} where {D, T}, Colon, UnitRange{Int}, AbstractVector{Int}}
@inline Base.getindex(d::DelaunayTriangulation, i::Int) = d.indices[i]
@inline Base.getindex(d::DelaunayTriangulation, i::tps) = d.indices[i]
@inline Base.getindex(d::DelaunayTriangulation, i::Int, j::tps) = d.indices[i, j]
@inline Base.getindex(d::DelaunayTriangulation, i::tps, j::tps) = d.indices[i, j]
@inline Base.getindex(d::DelaunayTriangulation, i::Int, j::Colon) = d.indices[i]
@inline Base.getindex(d::DelaunayTriangulation, i::tps, j::Colon) = d.indices[i]
@inline Base.getindex(d::DelaunayTriangulation, i::Colon, j::Int) = d.indices[i, j]
@inline Base.getindex(d::DelaunayTriangulation, i::Colon, j::Colon) = d.indices
@inline Base.getindex(d::DelaunayTriangulation, i::Colon, j::tps) = d.indices[i, j]
dimension(::DelaunayTriangulation{D,T}) where {D,T} = D
@inline Base.eltype(d::DelaunayTriangulation{D,T}) where {D,T} = T

#Base.unique(d::DelaunayTriangulation) = Base.unique(d.indices.data)
#Base.unique(d::DelaunayTriangulation, dims) = Base.unique(d.indices.data, dims)

#import Base: ==
#function ==(d₁::DelaunayTriangulation, d₂::DelaunayTriangulation)
#    d₁.indices == d₂.indices
#end

"""
    indices(d::DelaunayTriangulation, i::Int)

Return the indices corresponding to the i-th simplex of
the triangulation.
"""
indices(d::DelaunayTriangulation, i::Int) = d.indices[i]

"""
    indices(d::DelaunayTriangulation, i::Int, j::Int)

Return index corresponding to the j-th vertex of the
i-th simplex of the triangulation.
"""
indices(d::DelaunayTriangulation, i::Int, j::Int) = d.indices[i][j]


####################################
# Pretty printing.
####################################

function summarise(d::DelaunayTriangulation)
    _type = typeof(d)
    n_simplices = length(d)
    D = dimension(d)
    summary = "$_type with $n_simplices simplices\n"
    return join([summary, matstring(d.indices)], "")
end

Base.show(io::IO, r::DelaunayTriangulation) = println(io, summarise(r))


function delaunay(d::Dataset)
    triang = Dataset(delaunay(Matrix(d)))
    DelaunayTriangulation(triang)
end

function delaunay(r::Embedding)
    triang = delaunay(r.points)
    DelaunayTriangulation(triang)
end


export
delaunay,
DelaunayTriangulation,
indices

end
