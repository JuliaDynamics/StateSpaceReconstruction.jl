using Reexport
@reexport module Embeddings

using Parameters
using ..TimeSeries: SingleTimeSeries
using Simplices: Delaunay.delaunayn
using SimplexSplitting: centroids_radii2, heaviside0


abstract type Embedding end


"""
An embedding of a set of points.

`points::Array{Float64, 2}`
    The points furnishing the embedding

`ts::Vector{Vector{Float64}}`
    The time series used to construct the embedding. One for each column of `embedding`.

`ts_inds::Vector{Int}`
    Which time series are in which column of `embedding`?

`embedding_lags::Vector{Int}`
    Embedding lag for each column of `embedding`

`dim::Int`
    The dimension of the embedding
"""
@with_kw struct GenericEmbedding <: Embedding
    points::Array{Float64, 2} = Array{Float64, 2}(0, 0)
    ts::Vector{SingleTimeSeries{Float64}} = Vector{SingleTimeSeries{Float64}}(0)
    ts_inds::Vector{Int} = Int[]
    embedding_lags::Vector{Int} = Int[]
    dim::Int = 0
end

"""
An embedding in which the last point is guaranteed to lie within the convex
hull of the preceding points.

`points::Array{Float64, 2}`
    The points furnishing the embedding

`ts::Vector{Vector{Float64}}`
    The time series used to construct the embedding. One for each column of `embedding`.

`ts_inds::Vector{Int}`
    Which time series are in which column of `embedding`?

`embedding_lags::Vector{Int}`
    Embedding lag for each column of `embedding`

`dim::Int`
    The dimension of the embedding
"""
@with_kw struct InvariantEmbedding <: Embedding
    points::Array{Float64, 2} = Array{Float64, 2}(0, 0)
    ts::Vector{SingleTimeSeries{Float64}} = Vector{SingleTimeSeries{Float64}}(0)
    ts_inds::Vector{Int} = Int[]
    embedding_lags::Vector{Int} = Int[]
    dim::Int = 0
end



"""
	embed(ts::Vector{SingleTimeSeries{Float64}},
			ts_inds::Vector{Int},
			embedding_lags::Vector{Int})

Generalised embedding of SingleTimeSeries.
"""
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

    GenericEmbedding(
        points = E,
        ts = ts,
        ts_inds = ts_inds,
        embedding_lags = embedding_lags,
        dim = dim)
end

embed(ts::Vector{Vector{Float64}}, ts_inds::Vector{Int}, embedding_lags::Vector{Int}) =
    embed([SingleTimeSeries(ts[i]) for i = 1:length(ts)], ts_inds, embedding_lags)

embed(ts::Vector{Vector{Float64}}) = embed(
    	[SingleTimeSeries(ts[i]) for i = 1:length(ts)],
    	[i for i in 1:length(ts)],
    	[0 for i in 1:length(ts)])

embed(ts::Vector{Vector{Int}}, ts_inds::Vector{Int}, embedding_lags::Vector{Int}) =
    embed([SingleTimeSeries(float.(ts[i])) for i = 1:length(ts)], ts_inds, embedding_lags)

embed(ts::Vector{Vector{Int}}) = embed(
	[SingleTimeSeries(float.(ts[i])) for i = 1:length(ts)],
	[i for i in 1:length(ts)],
	[0 for i in 1:length(ts)])


include("embedding/invariantize.jl")


using RecipesBase

@recipe function f(E::GenericEmbedding)
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


@recipe function f(E::InvariantEmbedding)
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

export
    Embedding,
    GenericEmbedding, embed,
    InvariantEmbedding, invariantize, is_invariant_under_linearmap,
    plot

end
