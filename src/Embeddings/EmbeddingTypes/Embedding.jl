
include("EmbeddingData.jl")

"""
    Embedding{D, T}


## Fields
Assume `r` is a embedding, then

- **`r.points`**: The embedded state vectors, represented as an array of column
    vectors.
- **`r.embeddingdata`**: Information on what data and lags were used
    to construct the embedding.
"""
mutable struct Embedding{D, T} <: AbstractEmbedding{D, T}
    points::AbstractArray{T, 2}
    embeddingdata::EmbeddingData{D, T}
end

export Embedding
