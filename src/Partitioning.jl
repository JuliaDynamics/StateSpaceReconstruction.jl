using Reexport
@reexport module Partitioning

using Parameters
using ..Embeddings:
			Embedding,
			GenericEmbedding,
			InvariantEmbedding
using GroupSlices:
			groupslices,
			firstinds,
			groupinds
using Simplices: Delaunay.delaunayn
using SimplexSplitting:
			centroids_radii2,
			simplex_volumes,
			orientations


abstract type Partition end


include("partitioning/rectangularbinning.jl")
include("partitioning/triangulation.jl")


export
    Partition,
    Triangulation, triangulate,
    EquidistantBinning, bin_equidistant
end
