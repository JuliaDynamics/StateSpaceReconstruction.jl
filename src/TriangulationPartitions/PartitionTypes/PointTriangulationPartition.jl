
mutable struct PointTriangulationPartition{D, T} <: TriangulationPartition{D, T}
	data::AbstractArray{T, 2}
	simplexindices::DelaunayTriangulation
end


export PointTriangulationPartition
