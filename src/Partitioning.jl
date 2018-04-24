using Reexport
@reexport module Partitioning

using Parameters
using ..Embeddings:
			Embedding,
			GenericEmbedding,
			LinearlyInvariantEmbedding
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
	# Triangulation type and subtypes
	Triangulation,
    GenericTriangulation,
	LinearlyInvariantTriangulation,

	# Methods that dispatches on Triangulation subtypes
	triangulate,
	maybeintersecting_simplices,

	# Rectangular binning
    EquidistantBinning, bin_equidistant
end
