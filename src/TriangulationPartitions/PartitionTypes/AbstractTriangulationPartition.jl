
abstract type AbstractTriangulationPartition{D, T} end

Base.eltype(::AbstractTriangulationPartition{D,T}) where {D,T} = T

data(tp::AbstractTriangulationPartition) = tp.data
simplexindices(tp::AbstractTriangulationPartition) = tp.simplexindices
getsimplex(tp::AbstractTriangulationPartition, i) = tp.simplexindices[i]
nsimplices(tp::AbstractTriangulationPartition) = length(tp.simplexindices)
dimension(tp::AbstractTriangulationPartition) = dimension(tp.simplexindices)


# Indexing returns a tuple of vertices and indices.
#Base.index(tp::ATP, i) =

export
AbstractTriangulationPartition,
data,
simplexindices,
getsimplex,
nsimplices,
dimension
