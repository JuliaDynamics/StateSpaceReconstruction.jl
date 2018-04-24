"""
A triangulation of a cloud of embedded points into disjoint simplices.

The triangulation is mutable, so that the following will work:

```julia
# Triangulate a set of random points in 3D space.
t = triangulate(rand(20, 3))

# Refine triangulation until all simplices are below the mean radius of the original
# triangulation.
target_radius = mean(t.radii)
refine_variable_k!(t, target_radius)
```

"""

abstract type Triangulation <: Partition end

@with_kw mutable struct GenericTriangulation <: Triangulation
    embedding::Embedding = Embedding()
    # The vertices of the triangulation
    points::Array{Float64, 2} = Array{Float64, 2}(0, 0)

    # The image vertices of the triangulation
    impoints::Array{Float64, 2} = Array{Float64, 2}(0, 0)

    # Array of indices referencing the vertices furnishing each simplex, expressed both in terms
    # of the original points and their images under the linear forward map.
    simplex_inds::Array{Int, 2} = Array{Float64, 2}(0, 0)

    # Some properties of the simplices furnishing the triangulation
    centroids::Array{Float64, 2} = Array{Float64, 2}(0, 0)
    radii::Vector{Float64} = Float64[]
    centroids_im::Array{Float64, 2}  = Array{Float64, 2}(0, 0)
    radii_im::Vector{Float64} = Float64[]
    orientations::Vector{Float64} = Float64[]
    orientations_im::Vector{Float64} = Float64[]
    volumes::Vector{Float64} = Float64[]
    volumes_im::Vector{Float64} = Float64[]
end


@with_kw mutable struct LinearlyInvariantTriangulation <: Triangulation
    embedding::Embedding = Embedding()
    # The vertices of the triangulation
    points::Array{Float64, 2} = Array{Float64, 2}(0, 0)

    # The image vertices of the triangulation
    impoints::Array{Float64, 2} = Array{Float64, 2}(0, 0)

    #=
    Array of indices referencing the vertices furnishing each simplex,
    expressed both in terms of the original points and their images under the
    linear forward map.
    =#
    simplex_inds::Array{Int, 2} = Array{Float64, 2}(0, 0)

    # Some properties of the simplices furnishing the triangulation
    centroids::Array{Float64, 2} = Array{Float64, 2}(0, 0)
    radii::Vector{Float64} = Float64[]
    centroids_im::Array{Float64, 2}  = Array{Float64, 2}(0, 0)
    radii_im::Vector{Float64} = Float64[]
    orientations::Vector{Float64} = Float64[]
    orientations_im::Vector{Float64} = Float64[]
    volumes::Vector{Float64} = Float64[]
    volumes_im::Vector{Float64} = Float64[]
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

    LinearlyInvariantTriangulation(
        embedding = E,
        points = points,
        impoints = impoints,
        simplex_inds = simplex_inds,
        centroids = c,
        radii = r,
        volumes = vol,
        centroids_im = c_im,
        radii_im = r_im,
        volumes_im = vol_im,
        orientations = o,
        orientations_im = o_im)
end



function triangulate(E::GenericEmbedding)
    points = E.points[1:end-1, :]
    simplex_inds = delaunay_triang(points)
    impoints = E.points[2:end, :]
    c, r = centroids_radii2(points, simplex_inds)
    c_im, r_im = centroids_radii2(impoints, simplex_inds)
    vol = simplex_volumes(points, simplex_inds)
    vol_im = simplex_volumes(impoints, simplex_inds)
    o = orientations(points, simplex_inds)
    o_im = orientations(impoints, simplex_inds)

    GenericTriangulation(
        embedding = E,
        points = points,
        impoints = impoints,
        simplex_inds = simplex_inds,
        centroids = c,
        radii = r,
        volumes = vol,
        centroids_im = c_im,
        radii_im = r_im,
        volumes_im = vol_im,
        orientations = o,
        orientations_im = o_im)
end

triangulate(pts::AbstractArray{Float64, 2}) = triangulate(embed(pts))



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
