if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end


using StateSpaceReconstruction
using Base.Test

include("embedding_colwise.jl")
include("partition.jl")
