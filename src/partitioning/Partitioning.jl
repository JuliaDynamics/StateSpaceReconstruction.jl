using Reexport
@reexport module Partitioning

using Parameters
using Plots
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

abstract type Partition end

# This functionality should move to Simplices.jl in the future
include("../misc/simplexoperations.jl")

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
point_representatives,
unique_rows_info,

# Rectangular binning
AbstractRectangularBinning,
RectangularBinning,
dimension,
npoints,
assign_bin_labels,
assign_coordinate_labels,
minima_and_stepsizes,
plot_partition,
marginal_visitation_freq


end
