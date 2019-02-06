using LinearAlgebra
import Simplices.simplexintersection
import Simplices.intersectingvertices


################################
# Exports
################################
export
AbstractSimplex,

# Size
npoints,
nvertices,
dimension,

# Properties of the simplex
orientation,
volume,
centroid,
radius,
issingular,

# Generating new points/simplices from existing ones
generate_interior_points,
generate_interior_simplex,
generate_exterior_point,
generate_exterior_points,
generate_exterior_simplex,
generate_intersecting_simplex,

# Intersections between simplices
intersect, ∩,
intersect_vertices,

# Plotting
connectvertices,
splitaxes


################################################################################
# Abstract simplex types
################################################################################
abstract type AbstractSimplex{D, T} end
abstract type AbstractImmutableSimplex{D, T} <: AbstractSimplex{D, T} end
abstract type AbstractMutableSimplex{D, T} <: AbstractSimplex{D, T} end


#######################
# Indexing
#######################
# Return the i-th point as column vector
Base.getindex(s::AbstractSimplex, i) = s.vertices[i]
Base.getindex(s::AbstractSimplex, i, j) = s.vertices[i][j]
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

###########################################################
# Vector{vertex} representation => Array representation
###########################################################
Base.Array(s::AbstractSimplex) = hcat(s.vertices...,)

#######################
# Sizes
#######################
dimension(s::AbstractSimplex) = length(s[1])
npoints(s::AbstractSimplex) = length(s)
nvertices(s::AbstractSimplex) = length(s)

export dimension, npoints, nvertices

#########################
# Orientation and volumes
#########################
# Orientation (convention: append rows of ones at top)
function orientation(s::AbstractSimplex)
    det([ones(1, dimension(s) + 1); s[:, :]])
end

""" Get the unscaled volume of the simplex. Divide by factorial(dim) to get the true volume."""
volume(s::AbstractSimplex) = abs(orientation(s))

#######################
# Centroid
#######################
function centroid(s::AbstractSimplex)
    D = dimension(s)
    # Results in a dim-by-1 matrix, but we just want a vector, so drop the last dimension
    centroid = dropdims(s[:, :] * (ones(D + 1, 1) / (D + 1)), dims = 2)
end

#######################
# Radius
#######################
function radius(s::AbstractSimplex)
    D = dimension(s)
    centroidmatrix = repeat(centroid(s), 1, D + 1)
    radius = sqrt(maximum(ones(1, D) * ((s[:, :] .- centroidmatrix) .^ 2)))
end

"""
    issingular(simplex::AbstractArray{T, 2}) where {T<:Number}

Determines if a simplex is singular by checking if any of its vertices are
identical.
"""
function issingular(simplex::AbstractArray{T, 2}) where {T<:Number}
    length(unique(simplex.vertices)) != length(simplex.vertices)
end


############################################
# Generating interior points/simplices
############################################

"""
    generate_interior_points(polytope, n::Int)

Generate `n` points in the interior of the polytope `p`, which
is represented by a dim-by-(# points) array.
"""
function generate_interior_points(polytope, n::Int)
    R = rand(size(polytope, 2), n)
    convex_coeffs = ((1 ./ sum(R, dims = 1)) .* R)
    polytope * convex_coeffs
end

"""
    generate_interior_point(simplex::T) where {T <: AbstractSimplex}

Generate a point that lies in the interior of `simplex`.
"""
function generate_interior_point(simplex::T) where {T <: AbstractSimplex}
    S = Array(simplex)
    dim = dimension(simplex)
    R = rand(dim + 1)
    convex_coeffs = ((1 ./ sum(R, dims = 1)) .* R)
    S * convex_coeffs
end

"""
    generate_interior_points(simplex::T, n::Int) where {T <: AbstractSimplex}

Generate `n` points that lies in the interior of `simplex`.
"""
function generate_interior_points(simplex::T, n::Int) where {T <: AbstractSimplex}
    S = Array(simplex)
    dim = dimension(simplex)

    # Random linear combination coefficients
    R = rand(dim + 1, n)

    # Normalise the coefficients so that they sum to one. We can then create the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, dims = 1)) .* R
    inside_pts = S * normalised_coeffs
    [inside_pts[:, i] for i = 1:n]
end


"""
    generate_interior_simplex(simplex::T) where {T <: AbstractSimplex}

Generate a new simplex contained in the interior of `simplex`.
"""
function generate_interior_simplex(simplex::T) where {T <: AbstractSimplex}
    # Create n_pts random points contained within the simplex
    n_pts = length(simplex)
    interior_pts = generate_interior_points(simplex, n_pts)

    # Return a simplex of the same type as the input.
    T([interior_pts[:, i] for i = 1:n_pts])
end

"""
    generate_exterior_point(simplex::T) where {T <: AbstractSimplex}

Generate a point in the exterior of `simplex`.
"""
function generate_exterior_point(simplex::T) where {T <: AbstractSimplex}
    parentsimplex = Array(simplex)
    dim = dimension(simplex)

    # Random linear combination coefficients
    R = rand(dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, dims = 1)) .* R
    normalised_coeffs[1] += 1.5

    parentsimplex * normalised_coeffs
end

"""
    generate_exterior_points(simplex::T, n::Int) where {T <: AbstractSimplex}

Generate `n` points in the exterior of `simplex`.
"""
function generate_exterior_points(simplex::T, n::Int) where {T <: AbstractSimplex}
    [generate_exterior_point(simplex) for i = 1:n]
end

"""
    generate_exterior_simplex(simplex::T) where {T <: AbstractSimplex}

Generate a new simplex lying in the exterior of `simplex`.
"""
function generate_exterior_simplex(simplex::T) where {T <: AbstractSimplex}
    T(generate_exterior_points(simplex, dimension(simplex) + 1))
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

"""
    generate_intersecting_simplex(simplex::T) where {T <: AbstractSimplex}

Generate a new simplex that intersects with `simplex` in a nontrivial way.
"""
function generate_intersecting_simplex(simplex::T) where {T <: AbstractSimplex}
    dim = dimension(simplex)

    n_inside = rand(1:dim)
    n_outside = dim + 1 - n_inside

    verts_inside = generate_interior_points(simplex, n_inside)
    verts_outside = generate_exterior_points(simplex, n_outside)

    T([verts_inside; verts_outside])
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


################################
# Helper functions for plotting
################################
function connectvertices(s::T) where {T <: AbstractSimplex}
    # Convert to array
    s_arr = Array(s)
    hcat(s_arr, hcat([s_arr[:, i] for i in [1, 3, 1, 4, 2]]...,))
end

splitaxes(x) = ([x[k, :] for k = 1:size(x, 1)]...,)
