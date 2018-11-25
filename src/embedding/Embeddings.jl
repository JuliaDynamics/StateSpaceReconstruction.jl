using Reexport
@reexport module Embeddings

using StaticArrays
using Distributions
using Statistics

include("types.jl")
include("embed.jl")
include("delaunay.jl")
include("invariantize.jl")

end
