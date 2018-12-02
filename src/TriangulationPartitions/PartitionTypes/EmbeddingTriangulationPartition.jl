
mutable struct EmbeddingTriangulationPartition{D, T} <: TriangulationPartition{D, T}
	data::AbstractEmbedding
	simplexindices::DelaunayTriangulation
end

export EmbeddingTriangulationPartition
