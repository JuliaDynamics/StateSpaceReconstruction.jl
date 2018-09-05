
"""
    minima_and_stepsizes(points, ϵ) -> (Vector{Float}, Vector{Float})

Find the minima along each axis of the embedding, and computes appropriate
`stepsizes` given `ϵ`, which provide instructions on how to grid the space.
Assumes each point is a column vector.

Specifically, the binning procedure is controlled by the type of `ϵ`:

1. `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
2. `ϵ::Float` divides each axis into intervals of size `ϵ`.
3. `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
4. `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.
"""
function minima_and_stepsizes(points, ϵ)
    D = size(points, 1)
    n_pts = size(points, 2)

    axisminima = Vector{Float64}(D)
    top = Vector{Float64}(D)

    for i = 1:D
        axisminima[i] = minimum(points[i, :])
        top[i] = maximum(points[i, :])
    end
    axisminima = axisminima - (top - axisminima) / 100
    top = top + (top - axisminima) / 100

    stepsizes = Vector{Float64}(D)
    if typeof(ϵ) <: Float64
        stepsizes = [ϵ for i in 1:D]
    elseif typeof(ϵ) == Vector{Float64}
        stepsizes = ϵ
    elseif typeof(ϵ) <: Int
        stepsizes = (top - axisminima) / ϵ
    elseif typeof(ϵ) == Vector{Int}
        stepsizes = (top - axisminima) ./ ϵ
    end

    axisminima, stepsizes
end

"""
    assign_integer_bin_label_to_eachpoint!(
        A::Array{Int, 2},
        embedding_points::Array{Float64, 2},
        axisminima::Vector{Float64},
        stepsizes::Vector{Float64},
        npts::Int) -> Array{Int, 2}

Assign integer bin labels to the points of an embedding, storing the
labels in the preallocated array `A`.
"""
function assign_integer_bin_label_to_eachpoint!(
            A::Array{Int, 2},
            pts::Array{Float64, 2},
            axisminima::Vector{Float64},
            δs::Vector{Float64},
            npts::Int)
    @inbounds for i = 1:npts
        A[:, i] .= floor.(Int, abs(axisminima .- pts[:, i]) ./ δs)
    end
end

"""
    assign_bin_labels(points, ϵ)

Consider a rectangular grid specified by ϵ. Assign bin labels to the provided
`points` by checking which bins each point fall into. Each points is given a
label which is an integer encoding of the origin of the bin it falls into.

The following `ϵ` will work:

* `ϵ::Int` divide each axis into `ϵ` intervals of the same size.
* `ϵ::Float` divide each axis into intervals of size `ϵ`.
* `ϵ::Vector{Int}` divide the i-th axis into `ϵᵢ` intervals of the same size.
* `ϵ::Vector{Float64}` divide the i-th axis into intervals of size `ϵᵢ`.

The points are assumed to be provided as an array where each point is a column.
"""
function assign_bin_labels(points, ϵ)
    # Find the minima of the embedding and generate step
    # sizes along each axis
    bottom, stepsizes = minima_and_stepsizes(points, ϵ)
    if typeof(stepsizes) <: Union{Float64, Int}
        stepsizes = [stepsizes]
    end
    # How many points are there in the embedding, and what is its dimension?
    npts = size(points, 2)
    D = size(points, 1)

    # Each points of the embedding gets assigned to one bin.
    # The coordinates of each point of the original embedding are
    # assigned an integer number indicating which bin along the respective
    # dimension it falls into.
    visited_bins_inds = zeros(Int, D, npts)
    assign_integer_bin_label_to_eachpoint!(
        visited_bins_inds, points,
        bottom, stepsizes, npts)

    return visited_bins_inds
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
    assign_bin_labels(E.points, ϵ)
end
