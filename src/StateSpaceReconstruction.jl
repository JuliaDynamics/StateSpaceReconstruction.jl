module StateSpaceReconstruction

using Reexport
using Parameters

include("TimeSeries.jl")
include("Embeddings.jl")
include("Partitioning.jl")

# There might be many different types of state space reconstructions,
# depending on what type of analysis we want to run, so define a
# parent State Space Reconstruction (SSR) type
abstract type SSR end

struct GenericSSR
    embedding::Embedding
    partitions::Vector{Partition}
end



end # module
