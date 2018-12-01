using Reexport

@reexport module Triangulations
    include("Simplices/Simplices.jl")
    include("Delaunay/Delaunay.jl")
    include("PartitionTypes/TriangulationPartitionTypes.jl")
end
