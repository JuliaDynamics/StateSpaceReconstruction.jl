@reexport module RectangularPartitions

using Simplices: Delaunay
import ..Embeddings: AbstractEmbedding, Embedding
using ..GroupSlices: groupslices, firstinds, groupinds
using LinearAlgebra

abstract type RectangularPartition end

# This functionality should move to Simplices.jl in the future
include("bin_encoding.jl")
include("marginal_visitation_frequency.jl")

export
Partition,

# Rectangular binning
AbstractRectangularBinning,
RectangularBinning,
dimension,
npoints,
assign_bin_labels,
assign_coordinate_labels,
minima_and_stepsizes,
marginal_visitation_freq


end
