abstract type AbstractTriangulation end

abstract type AbstractTriangulationPartition{D, T} <: AbstractTriangulation end
abstract type AbstractMutableTriangulationPartition{D, T} <: AbstractTriangulationPartition{D, T} end

abstract type AbstractTriangulationPartitionFull{D, T} <: AbstractTriangulationPartition{D, T} end
abstract type AbstractMutableTriangulationPartitionFull{D, T} <: AbstractTriangulationPartitionFull{D, T} end

Base.eltype(::AbstractTriangulationPartition{D,T}) where {D,T} = T

####################
# Various
####################
data(tp::AbstractTriangulationPartition) = tp.data
simplexindices(tp::AbstractTriangulationPartition) = tp.simplexindices
getsimplex(tp::AbstractTriangulationPartition, i) = tp.simplexindices[i]
nsimplices(tp::AbstractTriangulationPartition) = length(tp.simplexindices)

dimension(tp::AbstractTriangulationPartition) = dimension(tp.simplexindices)


####################
# Pretty printing
####################
function summarise(tp::AbstractTriangulationPartition)
    _type = typeof(tp)
    n_simplices = nsimplices(tp)
    D = dimension(tp)
    summary = "$_type with $n_simplices simplices"
end

Base.show(io::IO, tp::AbstractTriangulationPartition) = println(io, summarise(tp))


####################
# Getting simplices
####################

"""
Materialize the simplices of a triangulation.
"""
function simplices(tpp::AbstractTriangulationPartition)

end

function getpoints(tpp::AbstractTriangulationPartition) end
function getsimplices(tpp::AbstractTriangulationPartition) end


export
AbstractTriangulationPartition,
AbstractMutableTriangulationPartition,
AbstractTriangulationPartitionFull,
AbstractMutableTriangulationPartitionFull,
data,
simplexindices,
getsimplex,
nsimplices,
dimension,
simplices,

getpoints,
getsimplices
