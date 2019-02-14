import DelayEmbeddings.Dataset

abstract type DatasetTriangulationPartition{D, T} <: AbstractTriangulationPartition{D, T} end
abstract type DatasetTriangulationPartitionFull{D, T} <: AbstractTriangulationPartitionFull{D, T} end

abstract type MutableDatasetTriangulationPartition{D, T} <: AbstractMutableTriangulationPartition{D, T} end
abstract type MutableDatasetTriangulationPartitionFull{D, T} <: AbstractMutableTriangulationPartitionFull{D, T} end

struct DatasetTriangulation{D, T} <: DatasetTriangulationPartition{D, T}
	data::Dataset
	simplexindices::DelaunayTriangulation
end

mutable struct MutableDatasetTriangulation{D, T} <: MutableDatasetTriangulationPartition{D, T}
	data::Dataset
	simplexindices::DelaunayTriangulation
end

# Full info
struct DatasetTriangulationFull{D, T} <: DatasetTriangulationPartitionFull{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
    simplices::Vector{Simplex}
    radii::Vector{Float64}
    orientations::Vector{Float64}
    volumes::Vector{Float64}
    centroids::Vector{Vector{Float64}}
end

# Full info
struct MutableDatasetTriangulationFull{D, T} <: MutableDatasetTriangulationPartitionFull{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
    simplices::Vector{Simplex}
    radii::Vector{Float64}
    orientations::Vector{Float64}
    volumes::Vector{Float64}
    centroids::Vector{Vector{Float64}}
end

PT = Union{DatasetTriangulationPartition, DatasetTriangulationPartitionFull,
    MutableDatasetTriangulationPartition, MutableDatasetTriangulationPartitionFull}
getpoints(pt::PT) = [pt.data[i] for i = 1:length(pt.data)]


getsimplices(pt::Union{DatasetTriangulationPartitionFull, MutableDatasetTriangulationPartitionFull}) =
    pt.simplices
getsimplices(pt::Union{DatasetTriangulationPartition, MutableDatasetTriangulationPartition}) =
    [Simplex(pt.data[pt.simplexindices[i]]) for i = 1:length(pt.simplexindices)]



####################
# Pretty printing
####################
function summarise(tp::Union{DatasetTriangulation, MutableDatasetTriangulation})
    _type = typeof(tp)
    n_simplices = nsimplices(tp)
    n_pts = length(tp.data)
    D = dimension(tp)
    summary = "$_type from $n_simplices simplices and $n_pts Datasets"
end

Base.show(io::IO, tp::Union{DatasetTriangulation, MutableDatasetTriangulation}) =
    println(io, summarise(tp))




export
DatasetTriangulation,
MutableDatasetTriangulation,
DatasetTriangulationFull,
MutableDatasetTriangulationFull
