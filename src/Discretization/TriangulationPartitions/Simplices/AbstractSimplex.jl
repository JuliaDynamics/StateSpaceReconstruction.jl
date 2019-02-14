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
intersect, 
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
# Subsampling simplices
#######################
include("simplex_subsampling.jl")

"""
    subsample(simplex::AbstractSimplex)

Subsample `simplex` with `n` points contained in its interior. If `sample_randomly = true`,
the points are randomly (uniformly) distributed within the simplex. 
If `sample_randomly = false` (default), then a shape-preserving simplex splitting 
routine is used to subdivide the simplex into subsimplices, and the centroids 
of these simplices are taken as the sampling points. For the latter approach,
the number of sampling points is *at least `n`*, but might be somewhat higher 
(the number of subsimplices generated by the splitting cannot be exactly `n`, 
so we use a splitting  scheme that yields a minimum of `n` subsimplices).
"""
function subsample(simplex::AbstractSimplex; n::Int = 100, sample_randomly::Bool = false)
    dim = length(simplex[1])

    # Convex expansion coefficients to generate the centroids of the subsimplices 
    # (if sample_randomly = false), or convex expansion coefficients to generate 
    # random points. 
    coeffs = subsample_coeffs(dim, n, sample_randomly)

    # Take the new points as as convex linear combinations of the vertices of 
    # the simplex
    pts = simplex[:, :] * coeffs
    [SVector{dim, Float64}(pts[:, i]) for i = 1:size(pts, 2)]
end

export subsample


"""
    subsample(simplex::AbstractSimplex, coeffs::Array{Int, 2})

Subsample `simplex` using the precomputed convex coefficients in `coeffs`, 
which contain instructions on how to combine the original vertices 
to create the subsampled points. Will return as many points 
as there are combinations/columns in `coeffs.`
"""
function subsample(simplex::AbstractSimplex, coeffs::Array{Float64, 2})
    dim = length(simplex[1])
    
    # Take the new points as as convex linear combinations of the vertices of 
    # the simplex
    pts = simplex[:, :] * coeffs
    [SVector{dim, Float64}(pts[:, i]) for i = 1:size(pts, 2)]
end

"""
    subsample(simplex1::AbstractSimplex, simplex2::AbstractSimplex;
        n::Int = 100, sample_randomly::Bool = false)

Subsample `simplex` with `n` points contained in its interior. 

If `sample_randomly = true`, the points are randomly (uniformly) distributed 
within the simplex. If `sample_randomly = false` (default), then a 
shape-preserving simplex splitting routine is used to subdivide the simplex 
into subsimplices, and the centroids 
of these simplices are taken as the sampling points. For the latter approach,
the number of sampling points is *at least `n`*, but might be somewhat higher 
(the number of subsimplices generated by the splitting cannot be exactly `n`, 
so we use a splitting  scheme that yields a minimum of `n` subsimplices).
"""
function subsample(simplex1::AbstractSimplex, simplex2::AbstractSimplex;
        n::Int = 100, sample_randomly::Bool = false)
    dim = length(simplex1[1])

    # Convex expansion coefficients to generate the centroids of the subsimplices 
    # (if sample_randomly = false), or convex expansion coefficients to generate 
    # random points. 
    coeffs = subsample_coeffs(dim, n, sample_randomly)

    # Take the new points as as convex linear combinations of the vertices of 
    # the simplex
    pts1 = simplex1[:, :] * coeffs
    pts2 = simplex2[:, :] * coeffs

    spts1 = [SVector{dim, Float64}(pts1[:, i]) for i = 1:size(pts1, 2)]
    spts2 = [SVector{dim, Float64}(pts2[:, i]) for i = 1:size(pts2, 2)]
    
    spts1, spts2
end

"""
    subsample(simplex::AbstractSimplex, coeffs::Array{Int, 2})

Subsample `simplex` using the precomputed convex coefficients in `coeffs`, 
which contain instructions on how to combine the original vertices 
to create the subsampled points. Will return as many points 
as there are combinations/columns in `coeffs.`
"""
function subsample(simplex1::AbstractSimplex, simplex2::AbstractSimplex, 
        coeffs::AbstractArray{Float64, 2})
    dim = length(simplex1[1])
    # Take the new points as as convex linear combinations of the vertices of 
    # the simplex
    pts1 = simplex1[:, :] * coeffs
    pts2 = simplex2[:, :] * coeffs

    spts1 = [SVector{dim, Float64}(pts1[:, i]) for i = 1:size(pts1, 2)]
    spts2 = [SVector{dim, Float64}(pts2[:, i]) for i = 1:size(pts2, 2)]
    
    spts1, spts2
end

#######################
# Refine simplices
#######################
include("simplex_splitting.jl")

function sumvertices(vertices, k::Int, dim::Int)
    s = zeros(Float64, dim) 
    for j = 1:k
        for i = 1:dim
            s[i] += vertices[j][i]
        end
    end
    return SVector{3, Float64}(s ./ k)
end


"""
    refine(simplex::Simplex; k::Int = 2)

Refine `simplex` with a shape-preserving sub-splitting of the simplex 
into ``k^{dim}`` subsimplices.
"""
function refine(simplex::AbstractSimplex; k::Int = 2,
        splitting_rules = nothing)

    dim = length(simplex.vertices) - 1
    
    # Rules for forming the strictly new vertices of the subtriangulation
    if splitting_rules == nothing
        rules, subtriangulation = simplex_splitting_rules(k, dim)
    elseif splitting_rules isa Tuple{Array{Int64,2},Array{Int64,2}}
        # all is good
        rules, subtriangulation = splitting_rules
    else
        throw(ArgumentError("Either specify both rules and subtriangulation, or neither of them."))
    end
    
    # How many new vertices are created each split?
    nverts_persplit = size(rules, 1)
    
    # The original vertices and strictly new sub-vertices formed by the splitting.
    all_vertices = Vector{SVector{dim, Float64}}(undef, nverts_persplit)
    
    # Generate the strictly new vertices for each subsimplex
    for j = 1:nverts_persplit
        # Pick the corresponding original vertices with indices contained in rules[j, :]
        # and compute their combination, which will be the j-th new vertex.
        all_vertices[j] = sumvertices(simplex.vertices[rules[j, :]], k, dim)
    end
    
    # Generate the new simplices
    subsimplices = Vector{Simplex{3, Float64}}(undef, k^dim)
    for i = 1:k^dim
        subsimplices[i] = Simplex(all_vertices[subtriangulation[:, i]])
    end
    
    return subsimplices
end



"""
    refine(simplex1::Simplex, simplex2::Simplex; k::Int = 2,
        splitting_rules::Union{Nothing, Tuple{Array{Float64, 2}, Array{Int, 2}}} = nothing)

Refine `simplex1` and `simplex2` with a shape-preserving sub-splitting of 
the simplices into ``k^{dim}`` subsimplices. 

This method may be used, for example, to keep track of which simplices 
correspond to each other when splitting a simplex and its image simplex 
(under the forward linear map of the vertices of the original simplex). 

`splitting_rules` may be generated beforehand to ensure that the same 
coefficients are used if many pairs of simplices are refined using 
the `simplex_splitting_rules` function.
"""
function refine(simplex1::AbstractSimplex, simplex2::AbstractSimplex; 
        k::Int = 2, 
        splitting_rules = nothing)

    dim = length(simplex1.vertices) - 1
    
    # Rules for forming the strictly new vertices of the subtriangulation
    if splitting_rules == nothing
        rules, subtriangulation = simplex_splitting_rules(k, dim)
    elseif splitting_rules isa Tuple{Array{Int64,2},Array{Int64,2}}
        # all is good
        rules, subtriangulation = splitting_rules
    else
        throw(ArgumentError("Either specify both rules and subtriangulation, or neither of them."))
    end
    
    # How many new vertices are created each split?
    nverts_persplit = size(rules, 1)
    
    # The original vertices and strictly new sub-vertices formed by the splitting.
    all_vertices1 = Vector{SVector{dim, Float64}}(undef, nverts_persplit)
    all_vertices2 = Vector{SVector{dim, Float64}}(undef, nverts_persplit)

    # Generate the strictly new vertices for each subsimplex
    for j = 1:nverts_persplit
        # Pick the corresponding original vertices with indices contained in rules[j, :]
        # and compute their combination, which will be the j-th new vertex.
        all_vertices1[j] = sumvertices(simplex1.vertices[rules[j, :]], k, dim)
        all_vertices2[j] = sumvertices(simplex2.vertices[rules[j, :]], k, dim)
    end
    
    # Generate the new simplices
    subsimplices1 = Vector{Simplex{3, Float64}}(undef, k^dim)
    subsimplices2 = Vector{Simplex{3, Float64}}(undef, k^dim)

    for i = 1:k^dim
        subsimplices1[i] = Simplex(all_vertices1[subtriangulation[:, i]])
        subsimplices2[i] = Simplex(all_vertices2[subtriangulation[:, i]])
    end
    
    return subsimplices1, subsimplices2
end

export refine 

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
