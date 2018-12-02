using Reexport
@reexport module Embeddings

import Simplices.Delaunay.delaunay

using StaticArrays
using Distributions
using Statistics
using NearestNeighbors


# Abstract embedding type
include("EmbeddingTypes/AbstractEmbedding.jl")

# Composite embedding types holding different types of data depending on
# the application.
include("EmbeddingTypes/Embedding.jl")
include("EmbeddingTypes/SimpleEmbedding.jl")
include("EmbeddingTypes/LinearlyInvariantEmbedding.jl")

# Invariantize embeddings (ensure last point lies within the convex hull
# of preceding points).
include("Invariantize/invariantize.jl")


# Interface with NearestNeighbors, allowing to create trees from embeddings.
include("NearestNeighbors/nearestneighbors.jl")

# Functions to perform delay embeddings
include("embed_functions/embed.jl")

end
