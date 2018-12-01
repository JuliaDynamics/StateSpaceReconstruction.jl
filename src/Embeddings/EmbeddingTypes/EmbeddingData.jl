
include("EmbeddingLags.jl")
include("EmbeddingPositions.jl")

"""
    EmbeddingData{N, T}

The data and parameters that goes into a embedding.

## Fields
If `rd` is a `EmbeddingData` instance, then

- **`rd.dataseries`**: gives a vector of dataseries (each itself being
a vector) that was used for the embedding.
- **`rd.in_which_pos`**: indicates which of the dataseries appear in which
columns of the embedding. For example, `rd.in_which_pos[2] = 3` means
that the 3rd timeseries (which can accessed by `rd.dataseries[3]`)
appears in the second column of the embedding.
- **`rd.at_what_lags[i]`**: gives the lag in the i-th column of the embedding.
- **`rd.labels[i]` returns the label of the i-th dataseries
- **`rd.descriptions[i]`**: returns the description of the i-th dataseries.
"""
struct EmbeddingData{D, T}
    dataseries::Vector{Union{AbstractVector{T}, SVector{D, T}}}
    labels::AbstractVector{String}
    descriptions::AbstractVector{String}
    in_which_pos::EmbeddingPositions
    at_what_lags::EmbeddingLags
end

####################################
# Indexing.
####################################
t = Union{Int, Colon, UnitRange, Vector{Int}}
Base.getindex(rl::EmbeddingLags, i::t) = Base.getindex(rl.lags, i)
Base.getindex(rp::EmbeddingPositions, i::t) = Base.getindex(rp.positions, i)

####################################
# Pretty printing.
####################################
function summarise(r::EmbeddingData)
    n  = length(r.dataseries)
    _type = typeof(r)
    """$_type consisting of $n data series."""
end

function description(r::EmbeddingData)
    S = String[]
    lags_s = string.(r.at_what_lags)
    positions_s = join(string.(r.in_which_pos))
    indexed_s = join(string.(1:length(r.dataseries)), " ")
    if all(r.labels .== "")

        push!(S, " The time series are both unnamed")
    else
        lbls = join([" The time series have labels", r.labels], ", ")
        push!(S, lbls)
    end
    push!(S, ", are indexed $indexed_s, and the data series appear")
    push!(S, " as columns of the embedding in $positions_s at $lags_s.")
    S = join([summarise(r); S], "")
    return S
end

Base.show(io::IO, r::EmbeddingData) = println(io, summarise(r))

export EmbeddingData
