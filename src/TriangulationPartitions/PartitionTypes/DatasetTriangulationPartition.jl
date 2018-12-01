import DynamicalSystemsBase.Dataset

mutable struct DatasetTriangulationPartition{D, T} <: TriangulationPartition{D, T}
	data::Dataset
	simplexindices::DelaunayTriangulation
end

export DatasetTriangulationPartition
