"""
   forwardlinearmap_invariant(pts::Array{T, 2}) where {T <: Number} -> Bool

Does the last point lie inside the convex hull of the preceding points?
"""
function forwardlinearmap_invariant(pts::Array{T, 2}) where {T <: Number}
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

    #=
    Find the row indices of the simplices that possibly contain the last point
    (meaning that dist(simplex_i, lastpoint) <= radius(simplex), so that the
    circumsphere of the simplex contains the last point.
    =#
    valid_simplex_indices = find(heaviside0(distdifferences) .*
                              collect(1:size(triangulation[2], 1)))

    n_validsimplices = size(valid_simplex_indices, 1)

    # Loop over valid simplices and check whether the corresponding simplex
    # contains the last point.
    i = 1
    lastpoint_contained = false

    #=
    Continue checking if simplices contain the last point until some simplex
    contains it. Then the point must necessarily be inside the convex hull of
    the triangulation formed by those simplices.
    =#
    while i <= n_validsimplices && !lastpoint_contained
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

        #=
        If the last convex expansion coefficient is positive, the last point is
        contained in the triangulation (because all the previous coefficients
        must have been nonnegative)
        =#
        if beta >= 0
            lastpoint_contained = true
        end

        i = i + 1
    end
    return lastpoint_contained
end

forwardlinearmap_invariant(E::T) where T<:StateSpaceReconstruction.AbstractEmbedding =
   forwardlinearmap_invariant(E.points)

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
function invariantize(emb::StateSpaceReconstruction.AbstractEmbedding;
                        max_increments = 20, verbose = false)
   pts = emb.points
   if size(unique(pts, dims = 1), 2) < size(pts, 2)
      @warn "Embedding points are not unique. Returning unmodified $emb."
      return emb
   end

   #=
   # Keep track of the embedding's center point and the original position of the
   # last point in the embedding, so we can move the last point along a line
   # from its original position towards the embedding's center, until the point
   # lies inside the convex hull of the preceding points.
   =#
   embedding_center = sum(pts, 1)/size(pts, 1)
   lastpoint = deepcopy(pts[end, :])
   dir = transpose(embedding_center) - lastpoint

   pts_removed = 0
   is_invariant = false
   percent_moved = 0.0
   dist_reduced = 0
   while !is_invariant && pts_removed <= max_increments
      if forwardlinearmap_invariant(pts)
         is_invariant = true
      else
         percent_moved += 1
         if verbose
            @warn """Moved last point $percent_moved % of the distance to convex hull origin."""
         end
         if isa(pts, Array)
            pts[end, :] = ceil(Int, lastpoint + dir*(percent_moved / 100))
         else
            pts[end, :] = lastpoint + dir*(percent_moved / 100)
         end
      end
   end

   if is_invariant
      return LinearlyInvariantEmbedding(
            pts[1:size(pts, 1), :],
            emb.which_ts,
            emb.in_which_pos,
            emb.at_what_lags,
            emb.dim
         )
   else
      @warn "Could not make embedding invariant. Returning empty embedding."
      return LinearlyInvariantEmbedding()
   end
end
