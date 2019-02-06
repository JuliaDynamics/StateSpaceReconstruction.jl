abstract type MatrixTriangulationPartition{D, T} <: AbstractTriangulationPartition{D, T} end
abstract type MatrixTriangulationPartitionFull{D, T} <: AbstractTriangulationPartitionFull{D, T} end

abstract type MutableMatrixTriangulationPartition{D, T} <: AbstractMutableTriangulationPartition{D, T} end
abstract type MutableMatrixTriangulationPartitionFull{D, T} <: AbstractMutableTriangulationPartitionFull{D, T} end

struct MatrixTriangulation{D, T} <: MatrixTriangulationPartition{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
end

mutable struct MutableMatrixTriangulation{D, T} <: MutableMatrixTriangulationPartition{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
end


# Full info
struct MatrixTriangulationFull{D, T} <: MatrixTriangulationPartitionFull{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
    simplices::Vector{Simplex}
    radii::Vector{Float64}
    orientations::Vector{Float64}
    volumes::Vector{Float64}
    centroids::Vector{Vector{Float64}}
end

# Full info
struct MutableMatrixTriangulationFull{D, T} <: MutableMatrixTriangulationPartitionFull{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
    simplices::Vector{Simplex}
    radii::Vector{Float64}
    orientations::Vector{Float64}
    volumes::Vector{Float64}
    centroids::Vector{Vector{Float64}}
end

MT = Union{MatrixTriangulationPartition, MatrixTriangulationPartitionFull,
    MutableMatrixTriangulationPartition, MutableMatrixTriangulationPartitionFull}
getpoints(mt::MT) = [mt.data[:, i] for i = 1:maximum(size(mt.data))]

getsimplices(mt::Union{MatrixTriangulationPartitionFull, MutableMatrixTriangulationPartitionFull}) =
    mt.simplices
getsimplices(mt::Union{MatrixTriangulationPartition, MutableMatrixTriangulationPartition}) =
    [Simplex(mt.data[:, mt.simplexindices[i]]) for i = 1:length(mt.simplexindices)]

####################
# Pretty printing
####################
function summarise(tp::Union{MatrixTriangulation, MutableMatrixTriangulation})
    _type = typeof(tp)
    n_simplices = nsimplices(tp)
    n_pts = maximum(size(tp.data))
    D = dimension(tp)
    summary = "$_type from $n_simplices simplices and $n_pts points\n"
end

Base.show(io::IO, tp::Union{MatrixTriangulation, MutableMatrixTriangulation}) =
    println(io, summarise(tp))



export
MatrixTriangulation,
MutableMatrixTriangulation,
MatrixTriangulationFull,
MutableMatrixTriangulationFull,
getpoints
