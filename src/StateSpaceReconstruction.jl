__precompile__(true)

module StateSpaceReconstruction

using Reexport
using StaticArrays
using LinearAlgebra

include("GroupSlices.jl")

include("Embeddings/Embeddings.jl")
include("RectangularPartitions/RectangularPartitions.jl")
include("TriangulationPartitions/TriangulationPartitions.jl")


# Rely on multiple dispatch for clashing function definitions in the submodules
dimension(E::AbstractEmbedding) = Embeddings.dimension(E)
npoints(E::AbstractEmbedding) = Embeddings.npoints(E)

dimension(s::Simplex) = Triangulations.Simplices.dimension(s)
npoints(s::Simplex) = Triangulations.Simplices.npoints(s)

dimension(dt::DelaunayTriangulation) = Triangulations.Delaunay.dimension(dt)

export dimension, npoints


end # module
