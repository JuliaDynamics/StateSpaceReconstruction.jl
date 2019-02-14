
using LinearAlgebra

import ..Embeddings: AbstractEmbedding, Embedding
using ..GroupSlices: groupslices, firstinds, groupinds

include("../plot_recipes/helperfunctions_rectgrid.jl")

include("bin_encoding.jl")
include("marginal_visitation_frequency.jl")


struct RectangularPartition{D, N}
    rectangles::Vector{Rectangle}
end

Base.length(r::RectangularPartition) = length(r.rectangles)
Base.getindex(r::RectangularPartition, i) = r.rectangles[i]
Base.firstindex(r::RectangularPartition, i) = 1
Base.lastindex(r::RectangularPartition, i) = Base.length(r.rectangles)
Base.iterate(r::RectangularPartition, state = 1) = iterate(r.rectangles, state)


export RectangularPartition