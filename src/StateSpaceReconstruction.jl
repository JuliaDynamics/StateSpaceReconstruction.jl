module StateSpaceReconstruction

using Reexport
using Parameters

include("TimeSeries.jl")
include("Embeddings.jl")
include("Partitioning.jl")


dimension(e::T where T<:AbstractEmbedding) = size(e.points, 2)
npoints(e::T where T<:AbstractEmbedding) = size(e.points, 1)
points(e::T where T <:AbstractEmbedding) = e.points
ntimeseries(e::T where T<:AbstractEmbedding) = length(e.which_ts)
timeseries(e::T where T<:AbstractEmbedding) = e.which_ts
which_ts(e::T where T<:AbstractEmbedding) = e.which_ts
in_which_pos(e::T where T<:AbstractEmbedding) = e.in_which_pos
at_what_lags(e::T where T<:AbstractEmbedding) = e.at_what_lags


# There might be many different types of state space reconstructions,
# depending on what type of analysis we want to run, so define a
# parent State Space Reconstruction (SSR) type
abstract type SSR end

struct GenericSSR
    embedding::Embedding
    partitions::Vector{Partition}
end



end # module
