__precompile__(true)

module StateSpaceReconstruction

using Reexport
using Parameters

include("embedding/Embeddings.jl")
include("partitioning/Partitioning.jl")

end # module
