using Reexport
@reexport module Embeddings

import Simplices.Delaunay.delaunay

using StaticArrays
using Distributions
using Statistics
using NearestNeighbors


# Abstract embedding type
include("AbstractEmbedding.jl")

# Composite embedding types holding different types of data depending on
# the application.
include("composite_types/Embedding.jl")
include("composite_types/SimpleEmbedding.jl")
include("composite_types/LinearlyInvariantEmbedding.jl")

# Invariantize embeddings (ensure last point lies within the convex hull
# of preceding points).
include("invariantizing/invariantize.jl")

# Interface with NearestNeighbors, allowing to create trees from embeddings.
include("NearestNeighbors/nearestneighbors.jl")

# Functions to perform delay embeddings
include("embed_functions/customembed.jl")

end
