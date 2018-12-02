using Reexport

@reexport module TriangulationPartitions
include("AbstractTriangulationPartition.jl")

using ..Delaunay
using ...Embeddings

abstract type TriangulationPartition{D, T} <: AbstractTriangulationPartition{D, T} end

include("DatasetTriangulationPartition.jl")
include("EmbeddingTriangulationPartition.jl")
include("PointTriangulationPartition.jl")


####################
# Pretty printing
####################
function summarise(tp::TriangulationPartition)
    _type = typeof(tp)
    n_simplices = nsimplices(tp)
    D = dimension(tp)
    #TODO: npts
    summary = "$D-dimensional $_type with $n_simplices simplices constructed from $n_pts points\n"
end

Base.show(io::IO, tp::TriangulationPartition) = println(io, summarise(tp))

export TriangulationPartition

end # module
