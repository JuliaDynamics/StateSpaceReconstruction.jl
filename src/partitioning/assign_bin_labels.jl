"""
    minima_and_stepsizes(E::AbstractEmbedding, ϵ) -> (Vector{Float}, Vector{Float})

Find the minima along each axis of the embedding,
and computes appropriate `stepsizes` given an `ϵ`.
"""
function minima_and_stepsizes(E::AbstractEmbedding, ϵ)
    emb = transpose(E.points)

    D = size(emb, 1) # dimension
    n_pts = size(emb, 2)

    bottom = [minimum(emb[:, i]) for i in 1:D]
    top = [maximum(emb[:, i]) for i in 1:D]
    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100

    stepsizes = Vector{Float64}(D)
    if typeof(ϵ) <: Float64
        stepsizes .= (top - bottom) * ϵ
    elseif typeof(ϵ) == Vector{Float64}
        stepsizes .= ϵ
    elseif typeof(ϵ) <: Int
        stepsizes .= (top - bottom) / ϵ
    elseif typeof(ϵ) == Vector{Int}
        stepsizes .= (top - bottom) ./ ϵ
    end

    bottom, stepsizes
end

"""
    assign_integer_bin_label_to_eachpoint!(
        A::Array{Int, 2},
        embedding_points::Array{Float64, 2},
        bottom::Vector{Float64},
        stepsizes::Vector{Float64},
        npts::Int) -> Array{Int, 2}

Assign integer bin labels to the points of an embedding, storing the
labels in the preallocated array `A`.
"""
function assign_integer_bin_label_to_eachpoint!(
            A::Array{Int, 2},
            pts::Array{Float64, 2},
            bottom::Vector{Float64},
            δs::Vector{Float64},
            npts::Int)
    @inbounds for i = 1:npts
        A[:, i] .= floor.(Int, abs.((pts[:, i] .- bottom)) ./ δs)
    end
end

"""
    assign_bin_labels(E::AbstractEmbedding, ϵ) -> Array{Int, 2}

Find bins visited by the reconstructed orbit and assign
to them integer tuples uniquely identifying the bins.

Each point of the embedding gets assigned one integer tuple label,
stored as column vectors in the returned array. If some bin
is visited more than once, there will be repeated columns.

# Details
Given a bin size ϵ, whics is Union{Int, Float} for cubic bins,
and Union{Vector{Int}, Vector{Float64}} for rectangular bins,
assign a tuple of integers uniquely identifying that bin.

Take, for example, the pair `(pᵢ, (n₁, n₂, n₃)`. Having a
3-element tuple associated with pᵢ means that the reconstruction
is 3D, so we have three coordinate axes - x₁, x₂ and x₃ - to
consider.

This generic example should be read as:
for point `pᵢ`, one must jump n₁ steps along x₁, n₂ steps
along x₂ and n₃ steps along x₃ (in units of ϵ₁, ϵ₂ and ϵ₃)
from the minima of each axis to reach the coordinates of `pᵢ`.

Each visited bin is thus uniquely identified by a tuple of
integers. `visited_bins_indices` assigns one such tuple
to each of the points in the embedding. These are gathered
in an array where column vectors represent the tuples
associated to each point.
"""
function assign_bin_labels(E::AbstractEmbedding, ϵ)
    pts = transpose(E.points)
    # Find the minima of the embedding and generate step
    # sizes along each axis
    bottom, stepsizes = minima_and_stepsizes(E, ϵ)

    # How many points are there in the embedding, and what is its dimension?
    npts = size(E.points, 1)
    D = E.dim

    # Each points of the embedding gets assigned to one bin.
    # The coordinates of each point of the original embedding are
    # assigned an integer number indicating which bin along the respective
    # dimension it falls into.
    visited_bins_inds = zeros(Int, D, npts)

    assign_integer_bin_label_to_eachpoint!(
        visited_bins_inds, pts,
        bottom, stepsizes, npts)

    return visited_bins_inds
end
