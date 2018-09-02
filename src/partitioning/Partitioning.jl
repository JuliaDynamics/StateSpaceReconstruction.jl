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
<<<<<<< HEAD
maybeintersecting_imsimplices,
=======
point_representatives,

>>>>>>> experimental_columnmajor_triangulation
unique_rows_info,

# Rectangular binning
AbstractRectangularBinning,
RectangularBinning,
dimension,
npoints,
assign_bin_labels,
marginal_visitation_freq


end
