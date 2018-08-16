using Reexport
@reexport module Partitioning

using Parameters
using ..Embeddings:
			AbstractEmbedding,
			Embedding,
			LinearlyInvariantEmbedding,
			SimpleEmbedding
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
	AbstractTriangulation,
    Triangulation,
	LinearlyInvariantTriangulation,

	# Methods that dispatches on Triangulation subtypes
	triangulate,
	maybeintersecting_simplices,

    unique_rows_info,

	# Rectangular binning
	RectangularBinning,
	bin_rectangular,
	    dimension, npoints, lowerbound, upperbound, stepsizes, indices_nonempty_bins,
	    unique_nonempty_bins, firstinds, groupinds, allinds

end
