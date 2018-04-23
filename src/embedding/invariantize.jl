"""
An embedding in which the last point is guaranteed to lie within the convex
hull of the preceding points.

`points::Array{Float64, 2}`
    The points furnishing the embedding

`ts::Vector{Vector{Float64}}`
    The time series used to construct the embedding. One for each column of `embedding`.

`ts_inds::Vector{Int}`
    Which time series are in which column of `embedding`?

`embedding_lags::Vector{Int}`
    Embedding lag for each column of `embedding`

`dim::Int`
    The dimension of the embedding
"""
@with_kw struct InvariantEmbedding <: Embedding
    points::Array{Float64, 2} = Array{Float64, 2}(0, 0)
    ts::Vector{SingleTimeSeries{Float64}} = Vector{SingleTimeSeries{Float64}}(0)
    ts_inds::Vector{Int} = Int[]
    embedding_lags::Vector{Int} = Int[]
    dim::Int = 0
end

"""
    is_invariant_under_linearmap(embedding::Embedding)

Does the last point of the embedding fall inside the convex hull of the
preceding points?
"""
function is_invariant_under_linearmap(pts::AbstractArray{Float64, 2})
    lastpoint = pts[end, :]
    dim = size(pts, 2)

    # Triangulate the embedding using all points but the last
    triangulation = delaunayn(pts[1:end-1, :])

    points = pts[1:end-1, :]
    simplex_indices = triangulation

    # Centroids and radii of simplices in the triangulation
    centroids, radii = centroids_radii2(points, simplex_indices)

    lastpoint_matrix = repmat(lastpoint', size(centroids, 1), 1)

    # Find simplices that can contain the last point (not all can)
    dists_lastpoint_and_centroids = sum((lastpoint_matrix - centroids).^2, 2)
    distdifferences = radii.^2 - dists_lastpoint_and_centroids

    # Find the row indices of the simplices that possibly contain the last point (meaning that
    # dist(simplex_i, lastpoint) <= radius(simplex), so the circumsphere of the simplex
    # contains the last point.
    valid_simplex_indices = find(heaviside0(distdifferences) .* collect(1:size(triangulation[2], 1)))

    n_validsimplices = size(valid_simplex_indices, 1)

    # Loop over valid simplices and check whether the corresponding simplex actually
    # contains the last point.
    i = 1
    lastpoint_contained = false

    # Continue checking if simplices contain the last point until some simplex contains it.
    # Then the point must necessarily be inside the convex hull of the triangulation formed
    # by those simplices.
    while i <= n_validsimplices && !lastpoint_contained
        # Subsample the valid simplex and compute its orientation
        simplex = pts[simplex_indices[valid_simplex_indices[i], :], :]
        orientation_simplex = det([ones(dim + 1, 1) simplex])

        beta = 1 # Convex expansion coefficient

        j = 1
        while j <= dim + 1 && beta >= 0
            tmp = copy(simplex)
            tmp[j, :] = copy(lastpoint)
            beta = det([ones(dim + 1, 1) tmp]) * sign(orientation_simplex)
            j = j + 1
        end

        # If the last convex expansion coefficient is positive, the last point is contained
        # in the triangulation (because all the previous coefficients must have been nonnegative)
        if beta >= 0
            lastpoint_contained = true
        end

        i = i + 1
    end
    # If the last point is contained in the triangulation, the set is invariant.
    return lastpoint_contained
end



"""
   invariantize_embedding(
      embedding::Array{Float64, 2};
      max_point_remove::Int = ceil(Int, size(embedding, 1)*0.05)
      )

If `remove_points = true`, iteratively remove the last point of an embedding
until it is invariant. If it is not possible to render the embedding invariant,
return an empty array. The default is to try to remove a maximum of ~5% of the
points of the original embedding before giving up. The number of points we're
allowed to remove can be set by providing the named argument `max_point_remove`.

If `remove_points = false`, incrementally move last point of the embedding
towards the origin until it lies within the convex hull of all preceding points.
"""
function invariantize(emb::GenericEmbedding; max_increments = 20)
   pts = emb.points
   if size(unique(pts, 1)) < size(pts)
      warn("Embedding points are not unique. Returning nothing.")
      return nothing
   end

   #=
   # Keep track of the embedding's centerpoint and the original position of the
   # last point in the embedding, so we can move the last point along a line
   # from its orinal position towards the embedding's center, until the point
   # lies inside the convex hull of the preceding points.
   =#
   embedding_center = sum(pts, 1)/size(pts, 1)
   lastpoint = pts[end, :]
   direction = embedding_center.' - lastpoint

   pts_removed = 0
   is_invariant = false
   percent_moved = 0.0

   while !is_invariant && pts_removed <= max_increments
      pts_removed += 1

      if is_invariant_under_linearmap(pts[1:(size(pts, 1) - pts_removed), :])
         is_invariant = true
      else
         percent_moved += 1
         warn("Moved point $percent_moved % towards origin to fit inside
               convex hull of previous points")
         pts[end, :] = lastpoint + direction*(percent_moved / 100)
      end
   end

   if is_invariant
      return InvariantEmbedding(
            points = pts[1:size(pts, 1), :],
            ts = emb.ts,
            ts_inds = emb.ts_inds,
            embedding_lags = emb.embedding_lags,
            dim = emb.dim
         )
   else
      warn("Could not make embedding invariant. Returning nothing.")
      return InvariantEmbedding()
   end
end
