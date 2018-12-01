
"""
     EmbeddingLags

Stores the information on which embedding lags were used to construct an
embedding.

## Fields
- **`lags`**: The lags which were used to construct the embedding. The i-th
    lag gives the lag for the i-th embedding variable.
"""
struct EmbeddingLags
    lags::Vector{Int}
end


export EmbeddingLags
