using Reexport

@reexport module Simplices
	import Distributions: Uniform
	
	# Abstract type and functions that should work on all composite types
	include("AbstractSimplex.jl")


	# Simplices whose vertices are represented by vectors
	include("composite_types/Simplex.jl")
	include("composite_types/MutableSimplex.jl")

	# Simplices whose vertices are represented by SVectors or MVectors
	include("composite_types/SSimplex.jl")
	include("composite_types/MutableSSimplex.jl")

	
	# Various functions
	include("point_contained_in_simplex.jl")

	# Plot recipes
	include("plotting/plot_recipes.jl")

end # module
