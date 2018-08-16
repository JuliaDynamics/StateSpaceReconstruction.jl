"""
A triangulation of a cloud of embedded points into disjoint simplices.

```julia
# Triangulate a set of random points in 3D space.
t = triangulate(rand(20, 3))

# Refine triangulation until all simplices are below the mean radius of the original
# triangulation.
target_radius = mean(t.radii)
refine_variable_k!(t, target_radius)
```

"""
abstract type AbstractTriangulation <: Partition end


#embedding(t::AbstractTriangulation) = t.embedding
#dimension(t::AbstractTriangulation) = dimension(t.embedding)
#npoints(t::AbstractTriangulation) = npoints(t.embedding)
#radii(t::AbstractTriangulation) = t.radii
#volumes(t::AbstractTriangulation) = t.volumes
#orientations(t::AbstractTriangulation) = t.orientations

"""
    `Triangulation <: AbstractTriangulation`

A triangulation type. Has the following fields:
1. `embedding::Embedding`.
2. `points::Array{Float64, 2}`. The vertices of the triangulation (a subset of the points
    of the embedding).
3. `impoints::Array{Float64, 2}`. The image vertices of the triangulation.
4. `simplex_inds::Array{Int, 2}`. # Array of indices referencing the vertices furnishing
    each simplex, expressed both in terms of the original points and their images under
    the linear forward map.
5. `centroids::Array{Float64, 2}`. Centroids of the simplices.
6. `radii::Vector{Float64}`. Radii of the simplices.
7. `centroids_im::Array{Float64, 2}`. Centroids of the image simplices.
8. `radii_im::Vector{Float64}`. Radii of the image simplices.
9. `orientations::Vector{Float64}`. Orientations of the simplices.
10. `orientations_im::Vector{Float64}`. Orientations of the image simplices.
11. `volumes::Vector{Float64}`. Volumes of the simplices.
12. `volumes_im::Vector{Float64}`. Volumes of the image simplices.
"""
struct Triangulation <: AbstractTriangulation
    embedding::Embedding
    points::Array{Float64, 2}
    impoints::Array{Float64, 2}
    simplex_inds::Array{Int, 2}
    centroids::Array{Float64, 2}
    centroids_im::Array{Float64, 2}
    radii::Vector{Float64}
    radii_im::Vector{Float64}
    orientations::Vector{Float64}
    orientations_im::Vector{Float64}
    volumes::Vector{Float64}
    volumes_im::Vector{Float64}
end



export dimension, npoints, radii, volumes, orientations

"""
    LinearlyInvariantTriangulation <: Triangulation`
A triangulation for which we have made sure the point corresponding to the last time
    index falls within the convex hull of the other points.
"""
struct LinearlyInvariantTriangulation <: AbstractTriangulation
    embedding::Embedding
    points::Array{Float64, 2}
    impoints::Array{Float64, 2}
    simplex_inds::Array{Int, 2}
    centroids::Array{Float64, 2}
    centroids_im::Array{Float64, 2}
    radii::Vector{Float64}
    radii_im::Vector{Float64}
    orientations::Vector{Float64}
    orientations_im::Vector{Float64}
    volumes::Vector{Float64}
    volumes_im::Vector{Float64}
end


"""
    triangulate(points::Array{Float64, 2})

Triangulate a set of vertices in N dimensions. `points` is an array of vertices,
where each row of the array is a point.
"""
function delaunay_triang(points::Array{Float64, 2})
    indices = delaunayn(points)
    return indices
end

function triangulate(E::LinearlyInvariantEmbedding)
    points = E.points[1:end-1, :]
    simplex_inds = delaunay_triang(points)
    impoints = E.points[2:end, :]
    c, r = centroids_radii2(points, simplex_inds)
    c_im, r_im = centroids_radii2(impoints, simplex_inds)
    vol = simplex_volumes(points, simplex_inds)
    vol_im = simplex_volumes(impoints, simplex_inds)
    o = orientations(points, simplex_inds)
    o_im = orientations(impoints, simplex_inds)

    LinearlyInvariantTriangulation(E,
        points,
        impoints,
        simplex_inds,
        c,
        c_im,
        r,
        r_im,
        o,
        o_im,
        vol,
        vol_im)
end

function triangulate(E::Embedding)
    points = E.points[1:end-1, :]
    simplex_inds = delaunay_triang(points)
    impoints = E.points[2:end, :]
    c, r = centroids_radii2(points, simplex_inds)
    c_im, r_im = centroids_radii2(impoints, simplex_inds)
    vol = simplex_volumes(points, simplex_inds)
    vol_im = simplex_volumes(impoints, simplex_inds)
    o = orientations(points, simplex_inds)
    o_im = orientations(impoints, simplex_inds)

    Triangulation(E,
        points,
        impoints,
        simplex_inds,
        c,
        c_im,
        r,
        r_im,
        o,
        o_im,
        vol,
        vol_im)
end

triangulate(pts::AbstractArray{Float64, 2}) = triangulate(embed(pts))

function Base.summary(t::T) where {T<:Triangulation}
    npts = size(t.embedding.points, 1)
    nsimplices = size(t.simplex_inds, 1)
    dim = size(t.embedding.points, 2)
    T = typeof(t.embedding)
    return "$dim-dimensional $(T) with $nsimplices simplices constructed from a $npts-pt $T"
end

function matstring(t::T) where {T<:Triangulation}
    fields = fieldnames(t)
    fields_str = String.(fields)
    maxlength = maximum([length(str) for str in fields_str]) + 2
    fields_str = [fields_str[i] *
                repeat(" ", maxlength - length(fields_str[i])) for i = 1:length(fields_str)]

    summaries = [join(":"*String(fields_str[i])*summary(getfield(t, fields[i]))*"\n") for i = 1:length(fields_str)] |> join
    infoline = "The following fields are available:\n"

    return summary(t)*"\n"*infoline*summaries
end

Base.show(io::IO, t::T) where {T<:Triangulation} = println(io, matstring(t))

"""
Find the indices of the simplices in the original triangulation that potentially
intersect with the image simplex with index `image_i`.
"""
function maybeintersecting_simplices(t::Triangulation, image_i::Int)
    inds_potential_simplices = Int[]

    n_simplices = length(t.radii)

    @inbounds for i = 1:n_simplices
        dist_difference = ((t.centroids_im[image_i] - t.centroids[i]).' *
                            (t.centroids_im[image_i] - t.centroids[i]) -
                                (t.radii_im[image_i] + t.radii[i])^2)[1]
        if dist_difference < 0
            push!(inds_potential_simplices, i)
        end
    end
    return inds_potential_simplices
end


"""
Find the indices of the image simplices in `t` that potentially intersect with
the original simplex with index `orig_i`.
"""
function maybeintersecting_imsimplices(t::Triangulation, orig_i::Int)
    inds_potential_simplices = Int[]

    n_simplices = length(t.radii)

    @inbounds for i = 1:n_simplices
        dist_difference = ((t.centroids[orig_i] - t.centroids_im[i]).' *
                            (t.centroids[orig_i] - t.centroids_im[i]) -
                                (t.radii[orig_i] + t.radii_im[i])^2)[1]
        if dist_difference < 0
            push!(inds_potential_simplices, i)
        end
    end
    return inds_potential_simplices
end
