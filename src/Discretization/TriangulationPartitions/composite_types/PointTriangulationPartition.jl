abstract type PointTriangulationPartition{D, T} <: AbstractTriangulationPartition{D, T} end
abstract type PointTriangulationPartitionFull{D, T} <: AbstractTriangulationPartitionFull{D, T} end

abstract type MutablePointTriangulationPartition{D, T} <: AbstractMutableTriangulationPartition{D, T} end
abstract type MutablePointTriangulationPartitionFull{D, T} <: AbstractMutableTriangulationPartitionFull{D, T} end

# Minimal
struct PointTriangulation{D, T} <: PointTriangulationPartition{D, T}
	data::Vector{Vector{T}}
	simplexindices::DelaunayTriangulation
end

mutable struct MutablePointTriangulation{D, T} <: MutablePointTriangulationPartition{D, T}
	data::Vector{Vector{T}}
	simplexindices::DelaunayTriangulation
end

# Full info
struct PointTriangulationFull{D, T} <: PointTriangulationPartitionFull{D, T}
	data::Vector{Vector{T}}
	simplexindices::DelaunayTriangulation
    simplices::Vector{Simplex}
    radii::Vector{Float64}
    orientations::Vector{Float64}
    volumes::Vector{Float64}
    centroids::Vector{Vector{Float64}}
end

# Full info
struct MutablePointTriangulationFull{D, T} <: MutablePointTriangulationPartitionFull{D, T}
	data::Vector{Vector{T}}
	simplexindices::DelaunayTriangulation
    simplices::Vector{Simplex}
    radii::Vector{Float64}
    orientations::Vector{Float64}
    volumes::Vector{Float64}
    centroids::Vector{Vector{Float64}}
end

PT = Union{PointTriangulationPartition, PointTriangulationPartitionFull,
    MutablePointTriangulationPartition, MutablePointTriangulationPartitionFull}
getpoints(pt::PT) = [pt.data[i] for i = 1:length(pt.data)]


getsimplices(pt::Union{PointTriangulationPartitionFull, MutablePointTriangulationPartitionFull}) =
    pt.simplices
getsimplices(pt::Union{PointTriangulationPartition, MutablePointTriangulationPartition}) =
    [Simplex(pt.data[pt.simplexindices[i]]) for i = 1:length(pt.simplexindices)]

####################
# Pretty printing
####################
function summarise(tp::Union{PointTriangulation, MutablePointTriangulation})
    _type = typeof(tp)
    n_simplices = nsimplices(tp)
    n_pts = length(tp.data)
    D = dimension(tp)
    summary = "$_type: $n_simplices simplices constructed from $n_pts points"
end

Base.show(io::IO, tp::Union{PointTriangulation, MutablePointTriangulation}) =
    println(io, summarise(tp))



export
PointTriangulation, MutablePointTriangulation,
PointTriangulationFull, MutablePointTriangulationFull
