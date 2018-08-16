using Reexport
@reexport module Embeddings

using RecipesBase
using Parameters

using Simplices: Delaunay.delaunayn
using SimplexSplitting: centroids_radii2, heaviside0

using ..TimeSeries: SingleTimeSeries

export
    AbstractEmbedding,
    Embedding,
    SimpleEmbedding,
    LinearlyInvariantEmbedding,
    embed,
    invariantize,
    is_invariant_under_linearmap

""" Abstract Embedding type. """
abstract type AbstractEmbedding end



"""
An embedding of a set of points. Has the fields
1. `points::Array{Float64, 2}`. The points furnishing the embedding
2. `ts::Vector{Vector{Float64}}`. The time series used to construct the embedding. One for
        each column of `embedding`.
3. `ts_inds::Vector{Int}`. Which time series are in which column of `embedding`?
4. `embedding_lags::Vector{Int}`. Embedding lag for each column of `embedding`
5. `dim::Int`. The dimension of the embedding
"""
struct Embedding <: AbstractEmbedding
    points::Array{Float64, 2}
    ts::Vector{SingleTimeSeries{Float64}}
    ts_inds::Vector{Int}
    embedding_lags::Vector{Int}
    dim::Int
end

"""
An embedding in which the last point is guaranteed to lie within the convex
hull of the preceding points.
"""
struct LinearlyInvariantEmbedding <: AbstractEmbedding
    points::Array{Float64, 2}
    ts::Vector{SingleTimeSeries{Float64}}
    ts_inds::Vector{Int}
    embedding_lags::Vector{Int}
    dim::Int
end

""" An embedding holding only its points and no information about the embedding itself."""
struct SimpleEmbedding <: AbstractEmbedding
    points::Array{Float64, 2}
end

dimension(e::T where T<:AbstractEmbedding) = size(e.points, 2)
npoints(e::T where T<:AbstractEmbedding) = size(e.points, 1)
points(e::T where T <:AbstractEmbedding) = e.points
ntimeseries(e::T where T<:AbstractEmbedding) = length(e.ts)
timeseries(e::T where T<:AbstractEmbedding) = e.ts
which_ts(e::T where T<:AbstractEmbedding) = e.ts
in_which_pos(e::T where T<:AbstractEmbedding) = e.ts_inds
at_what_lags(e::T where T<:AbstractEmbedding) = e.embedding_lags

export dimension, npoints, ntimeseries, timeseries, which_ts, in_which_pos, at_what_lags



function embed(ts::Vector{SingleTimeSeries{Float64}},
               ts_inds::Vector{Int},
               embedding_lags::Vector{Int})
    dim = length(ts_inds)
    minlag, maxlag = minimum(embedding_lags), maximum(embedding_lags)
    npts = length(ts[1].ts) - (maxlag + abs(minlag))
    E = zeros(Float64, npts, dim)

    for i in 1:length(ts_inds)
        ts_ind = ts_inds[i]
        TS = ts[ts_ind].ts
        lag = embedding_lags[i]

        if lag > 0
            E[:, i] = TS[((1 + abs(minlag)) + lag):(end - maxlag) + lag]
        elseif lag < 0
            E[:, i] = TS[((1 + abs(minlag)) - abs(lag)):(end - maxlag - abs(lag))]
        elseif lag == 0
            E[:, i] = TS[(1 + abs(minlag)):(end - maxlag)]
        end
    end

    Embedding(E, ts, ts_inds, embedding_lags, dim)
end

embed(ts::Vector{Vector{T}} where T<:Number) = embed(
    	[SingleTimeSeries(float.(ts[i])) for i = 1:length(ts)],
    	[i for i in 1:length(ts)],
    	[0 for i in 1:length(ts)]
)

"""
	embed(ts::Vector{SingleTimeSeries{Float64}},
			ts_inds::Vector{Int},
			embedding_lags::Vector{Int})

Embed a set of vectors.

## Arguments
1. `which_ts::Vector{Vector{Float64}}`. This is a vector containing the time series to embed.
    - Example: which_ts = [ts1, ts2].
2. `in_which_pos::Vector{Int}``. The length of in_which_pos gives the dimension of the
    embedding. The value of the ith element of in_which_pos indicates which time series in
    the ith column of the embedding.
    - **Example 1**: if `which_ts = [ts1, ts2]`, then we index ts1 as 1 and ts2 as 2.
        Setting `in_which_pos = [2, 2, 1]` will result in a 3-dimensional embedding where
        `ts2` will appear in columns 1 and 2, while `ts1` will appear in column 3.
    - **Example 2**: If `which_ts = [ts1, ts2, ts3]`, then `in_which_pos = [2,1,2,3,3]`
        results in a 5-dimensional embedding where `ts1`appears in column 2, `ts2` appears
        in columns 1 and 3, while `ts3`appears in columns 4 and 5.
3. `at_what_lags::Vector{Int}` sets the lag in each column. Must be the same length as
    `which_ts`.
    - **Example**: if `in_which_pos = [2, 2, 1]`, then  `at_what_lags = [1, 0, -1]` means
        that the lag in column 1 is 1, the lag in the second column is 0 and the lag in
        the third column is -1.

"""
embed(ts::Vector{Vector{T}} where T<:Number, ts_inds::Vector{Int}, embedding_lags::Vector{Int} where T<:Real) =
    embed([SingleTimeSeries(float.(ts[i])) for i = 1:length(ts)], ts_inds, embedding_lags)

"""
Default embedding of a `npts`-by-`dim` array of points.
"""
embed(A::AbstractArray{Float64, 2}) = embed(
    	[A[:, i] for i = 1:size(A, 2)],
    	[i for i in 1:size(A, 2)],
    	[0 for i in 1:size(A, 2)])

embed(A::AbstractArray{Int, 2}) = embed(float.(A))

embed(A::AbstractArray{Float64, 2}, ts_inds::Vector{Int}, embedding_lags::Vector{Int}) =
    embed([A[:, i] for i = 1:size(A, 2)], ts_inds, embedding_lags)

embed(A::AbstractArray{Int, 2}, ts_inds::Vector{Int}, embedding_lags::Vector{Int}) =
        embed([float.(A[:, i]) for i = 1:size(A, 2)], ts_inds, embedding_lags)

export embed

include("embedding/invariantize.jl")



@recipe function f(E::Embedding)
    if E.dim > 3
        warn("Embedding dim > 3, plotting three first axes")
        pts = E.points[:, 1:3]
    end
    pts = E.points
    X = pts[:, 1]
    Y = pts[:, 2]
    Z = pts[:, 3]
    X, Y, Z
end


end
