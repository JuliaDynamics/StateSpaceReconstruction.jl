using StateSpaceReconstruction
using Statistics: std
using LinearAlgebra
using Test

if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

include("embedding_colwise.jl")
include("partition.jl")
