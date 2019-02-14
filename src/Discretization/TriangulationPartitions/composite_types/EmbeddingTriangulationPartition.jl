import ..Embeddings: AbstractEmbedding

abstract type EmbeddingTriangulationPartition{D, T} <: AbstractTriangulationPartition{D, T} end
abstract type MutableEmbeddingTriangulationPartition{D, T} <: AbstractMutableTriangulationPartition{D, T} end

abstract type EmbeddingTriangulationPartitionFull{D, T} <: AbstractTriangulationPartitionFull{D, T} end
abstract type MutableEmbeddingTriangulationPartitionFull{D, T} <: AbstractMutableTriangulationPartitionFull{D, T} end

struct EmbeddingTriangulation{D, T} <: EmbeddingTriangulationPartition{D, T}
	data::AbstractEmbedding
	simplexindices::DelaunayTriangulation
end

mutable struct MutableEmbeddingTriangulation{D, T} <: MutableEmbeddingTriangulationPartition{D, T}
	data::AbstractEmbedding
	simplexindices::DelaunayTriangulation
end

# Full info
struct EmbeddingTriangulationFull{D, T} <: EmbeddingTriangulationPartitionFull{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
    simplices::Vector{Simplex}
    radii::Vector{Float64}
    orientations::Vector{Float64}
    volumes::Vector{Float64}
    centroids::Vector{Vector{Float64}}
end

# Full info
struct MutableEmbeddingTriangulationFull{D, T} <: MutableEmbeddingTriangulationPartitionFull{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
    simplices::Vector{Simplex}
    radii::Vector{Float64}
    orientations::Vector{Float64}
    volumes::Vector{Float64}
    centroids::Vector{Vector{Float64}}
end



MT = Union{EmbeddingTriangulationPartition, EmbeddingTriangulationPartitionFull,
    MutableEmbeddingTriangulationPartition, MutableEmbeddingTriangulationPartitionFull}
getpoints(mt::MT) = [mt.data[:, i] for i = 1:maximum(size(mt.data))]

getsimplices(mt::Union{EmbeddingTriangulationPartitionFull, MutableEmbeddingTriangulationPartitionFull}) =
    mt.simplices
getsimplices(mt::Union{EmbeddingTriangulationPartition, MutableEmbeddingTriangulationPartition}) =
    [Simplex(mt.data[:, mt.simplexindices[i]]) for i = 1:length(mt.simplexindices)]


####################
# Pretty printing
####################
function summarise(tp::Union{EmbeddingTriangulation, MutableEmbeddingTriangulation})
    _type = typeof(tp)
    n_simplices = nsimplices(tp)
    n_pts = size(tp.data, 2)
    D = dimension(tp)
    summary = "$_type from $n_simplices simplices and $n_pts points"
end

Base.show(io::IO, tp::Union{EmbeddingTriangulation, MutableEmbeddingTriangulation}) =
    println(io, summarise(tp))


export
EmbeddingTriangulation,
MutableEmbeddingTriangulation,
EmbeddingTriangulationFull,
MutableEmbeddingTriangulationFull
