__precompile__(true)

module StateSpaceReconstruction

using Reexport
using StaticArrays
using LinearAlgebra

include("GroupSlices.jl")

include("Embeddings/Embeddings.jl")
include("Discretization/RectangularPartitions/RectangularPartitions.jl")
include("Discretization/TriangulationPartitions/TriangulationPartitions.jl")


# Rely on multiple dispatch for clashing function definitions in the submodules
dimension(E::AbstractEmbedding) = Embeddings.dimension(E)
dimension(s::T) where {T<:AbstractSimplex} = TriangulationPartitions.Simplices.dimension(s)
dimension(dt::T) where {T<:AbstractDelaunayTriangulation} =
    TriangulationPartitions.Delaunay.dimension(dt)
npoints(s::T) where {T<:AbstractSimplex} = TriangulationPartitions.Simplices.npoints(s)
npoints(E::AbstractEmbedding) = Embeddings.npoints(E)


export dimension, npoints

end # module



# #
# #
# # function estimator_exact(pts::Vector{Vector{T}}) where T
# #
# #     DT = StateSpaceReconstruction.DelaunayTriangulation(pts)
# #
# #     # Tolerance for zero volume (if one of the simplices have nearly zero
# #     # volume, skip computing the intersection volume and just set it to
# #     # zero instead).
# #     n_simplices = length(DT.indices)
# #     ϵ::Float64 = 1e-5/n_simplices
# #
# #     dim = length(DT[1]) - 1
# #
# #     image_simplex = StateSpaceReconstruction.MutableSimplex([zeros(Float64, dim) for i = 1:(dim+1)])
# #     simplex = StateSpaceReconstruction.MutableSimplex([zeros(Float64, dim) for i = 1:(dim+1)])
# #
# #     transfermatrix = zeros(Float64, n_simplices, n_simplices) # intersecting volumes
# #
# #     println("Building transfer matrix")
# #     for i in 1:n_simplices
# #         [image_simplex[k] = pts[DT[i][k]] for k = 1:(dim + 1)]
# #         for j in 1:n_simplices
# #             [simplex[k] = pts[DT[j][k]] for j = 1:(dim + 1)]
# #
#             @show image_simplex
#         end
#     end
#
# end
#
#
#
# @time estimator_exact([rand(3) for i = 1:10])
#
#
#
# TO = zeros(Float64, n_simplices, n_simplices) # intersecting volumes
#
# for i in 1:n_simplices
#     imvol = t.volumes_im[i]
#     for j in 1:n_simplices
#         vol = t.volumes[j]
#         if vol * imvol > 0 && (vol/imvol) > ϵ
#             # Intersecting volume between these two simplices
#             TO[i, j] = simplexintersection(
#                             transpose(t.points[t.simplex_inds[j, :], :]),
#                             transpose(t.impoints[t.simplex_inds[i, :], :])
#                         ) / imvol
#         end
#     end
