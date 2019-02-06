if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end


using StateSpaceReconstruction
using Statistics: std
using LinearAlgebra
using Test

if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

include("Embeddings/test_Embeddings.jl")
include("Embeddings/test_invariantize.jl")

include("RectangularPartitions/test_RectangularPartitions.jl")
include("Triangulations/test_Simplices.jl")
include("Triangulations/test_Delaunay.jl")
include("Triangulations/test_Triangulations.jl")
