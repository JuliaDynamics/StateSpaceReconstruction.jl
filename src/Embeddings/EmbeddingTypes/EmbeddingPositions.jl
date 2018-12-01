
"""
     EmbeddingPositions

If the embedding was constructed from multivariate data, which of the
original data variables was mapped to which variable in the embedding?

## Fields
- **`positions`**: position[i] returns the index of the variable that goes
    into the i-th coordinate axis of the embedding.
"""
struct EmbeddingPositions
    positions::AbstractVector{Int}
end


export EmbeddingsPositions
