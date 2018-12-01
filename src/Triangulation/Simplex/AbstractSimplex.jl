using LinearAlgebra
import Simplices.simplexintersection
import Simplices.intersectingvertices

abstract type AbstractSimplex end



#######################
# Indexing
#######################
# Return the i-th point as column vector
Base.getindex(s::AbstractSimplex, i) = s.vertices[i]
Base.getindex(s::AbstractSimplex, i::Colon, j) = hcat(s.vertices[j]...,)
Base.getindex(s::AbstractSimplex, i::Colon, j::Colon) = hcat(s.vertices...,)

Base.length(s::AbstractSimplex) = length(s.vertices)
Base.size(s::AbstractSimplex) = length(s[1]), length(s)
Base.size(s::AbstractSimplex, i) = size(s)[i]
Base.IteratorSize(s::AbstractSimplex) = Base.HasLength()

Base.firstindex(s::AbstractSimplex) = 1
Base.lastindex(s::AbstractSimplex) = length(s)
Base.eachindex(s::AbstractSimplex) = Base.OneTo(length(s))
Base.iterate(s::AbstractSimplex, state = 1) = iterate(s.vertices, state)

#######################
# Sizes
#######################
dimension(s::AbstractSimplex) = length(s[1])
npoints(s::AbstractSimplex) = length(s)
nvertices(s::AbstractSimplex) = length(s)

#########################
# Orientation and volumes
#########################
# Orientation (convention: append rows of ones at top)
function orientation(s::AbstractSimplex)
    det([ones(1, dimension(s) + 1); s[:, :]])
end

volume(s::AbstractSimplex) = abs(orientation(s))/factorial(dimension(s))

#######################
# Centroid
#######################
function centroid(s::AbstractSimplex)
    D = dimension(s)
    centroid = s[:, :] * (ones(D + 1, 1) / (D + 1))
end

#######################
# Radius
#######################
function radius(s::AbstractSimplex)
    D = dimension(s)
    centroidmatrix = repeat(centroid(s), 1, D + 1)
    radius = sqrt(maximum(ones(1, D) * ((s[:, :] .- centroidmatrix) .^ 2)))
end


#######################
# Intersections
#######################
"""
    intersect(s1::AbstractSimplex, s2::AbstractSimplex)

Compute the volume intersection between two simplices.
"""
function Base.intersect(s1::AbstractSimplex, s2::AbstractSimplex)
    simplexintersection(s1[:, :], s2[:, :])
end

"""
    ∩(s1::AbstractSimplex, s2::AbstractSimplex)

Compute the volume intersection between two simplices.
"""
∩(s1::AbstractSimplex, s2::AbstractSimplex) = intersect(s1, s2)


"""
    intersect_vertices(s1::AbstractSimplex, s2::AbstractSimplex)

Return the vertices forming the convex hull of the volume intersection
between two simplices.
"""
function intersect_vertices(s1::AbstractSimplex, s2::AbstractSimplex)
    intersectingvertices(s1[:, :], s2[:, :])
end


#######################
# Pretty printing
#######################
function summarise(s::AbstractSimplex)
    _type = typeof(s)
    n_vertices = npoints(s)
    D = dimension(s)

    matrixrepresentation = sprint(io -> show(IOContext(io, :limit=>true),
                MIME"text/plain"(), s[:, :]))
    matrixrepresentation = join(split(matrixrepresentation, '\n')[2:end], '\n')

    summary = join(["$D-dimensional $_type with $n_vertices vertices", "\n",
                    matrixrepresentation])
end

Base.show(io::IO, s::AbstractSimplex) = println(io, summarise(s))


export
AbstractSimplex,
dimension,
npoints,
nvertices,
orientation,
volume,
centroid,
radius,
intersect,
∩,
intersect_vertices
