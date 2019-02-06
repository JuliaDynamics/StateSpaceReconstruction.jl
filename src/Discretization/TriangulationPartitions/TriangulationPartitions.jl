using Reexport

@reexport module TriangulationPartitions
    # Necessary modules
    include("Simplices/Simplices.jl")
    include("DelaunayTriangulations/DelaunayTriangulations.jl")

	dimension(DT::DelaunayTriangulation) = DelaunayTriangulations.dimension(DT)

    # Abstract types and functions that work on all subtypes
    include("AbstractTriangulationPartition.jl")

    # Composite types for different types of data
    include("composite_types/DatasetTriangulationPartition.jl")
    include("composite_types/EmbeddingTriangulationPartition.jl")
    include("composite_types/PointTriangulationPartition.jl")
	include("composite_types/MatrixTriangulationPartition.jl")

    # Functions to create triangulations from data
    include("triangulate.jl")

    # Plot triangulations
    include("plot_recipes/plot_recipes.jl")

end


"""
    TriangulationPartitions

A module handling partitioned formed by triangulations.
"""
TriangulationPartitions
