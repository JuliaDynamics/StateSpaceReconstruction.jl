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
using Simplices:
	Delaunay.delaunayn,
	childpoint
using SimplexSplitting:
			centroids_radii2,
			simplex_volumes,
			orientations


abstract type Partition end

include("rectangularbinning.jl")
include("triangulation.jl")
include("visitation_frequency_marginals.jl")

export
Partition,
# Triangulation type and subtypes
AbstractTriangulation,
Triangulation,
LinearlyInvariantTriangulation,

# Methods that dispatches on Triangulation subtypes
triangulate,
maybeintersecting_simplices,
maybeintersecting_imsimplices,
unique_rows_info,

# Rectangular binning
AbstractRectangularBinning,
RectangularBinning,
dimension,
npoints,
assign_bin_labels,
marginal_visitation_freq


end
