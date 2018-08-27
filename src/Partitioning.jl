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
include("partitioning/assign_bin_labels.jl")

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
AbstractRectangularBinning,
RectangularBinning,
assign_bin_labels,
bin_rectangular,
dimension,
npoints


end
