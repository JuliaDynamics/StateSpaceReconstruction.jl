"""
    visitation_freq(
        along_which_axes::Union{Vector{Int}, AbstractUnitRange{Int}},
        visited_bins_inds::Array{Int, 2},
        npts::Int) -> Vector{Float64}

Given an array of integer tuples (represented as column vectors of the array),
estimate marginal or joint visitation frequencies. The argument
`along_which_axes` controls which axes to take the marginal visitation
frequencies along.

Setting `along_which_axes` to a range 1:D, where
D is the dimension of the corresponding state space, corresponds to
taking the joint visitation frequency.

`npts` is the number of points
in the embedding, which is used as a normalization factor.
"""
function marginal_visitation_freq(
        along_which_axes::Union{Vector{Int}, AbstractUnitRange{Int}},
        visited_bins_inds::Array{Int, 2},
        npts::Int)

    # Assuming there are M unique bins (or, equivalently, integer tuples),
    # groupslices(visited_bins_inds, 2) returns a vector of indices,
    # in which each visited bin label is associated with the i-th unique
    # bin. For example, if slices = [1, 2, 3, 1, 2, 3, 7, ...], then points
    # p₁ and p₄ are repeated (have identical labels), p₂ and p₅ are repeated,
    # and p₃ and p₆ are repeated.
    #
    # Subsetting some of the rows are equivalent to estimating the marginal
    # along the corresponding axes. If along_which_axes = 1:n_axes, then we're
    # computing the joint frequency distribution.
    slices = groupslices(visited_bins_inds[along_which_axes, :], 2)

    # The groupinds function takes these indices and groups the indices
    # corresponding to repeated points. For the example in the comment
    # above, group_repeated_inds = [[1, 4], [2, 5], [3, 6],
    # [7, indices of potential points shared with p₇], ....].
    group_repeated_inds = groupinds(slices)

    # Computing the joint probability for pᵢ is now just a matter of counting
    # the number of times each visited bin is visited. The number of times
    # the ith visited bin is visited by the orbit is then simply counting
    # how many elements there are in group_repeated_inds[i].
    n_visited_states = length(group_repeated_inds)

    # We'll now loop over one bit at a time, counting how many times
    # the bin is visited, then get visitation frequency for that bin
    # by dividing the number of times it is visited by the total
    # number of points.
    m = Vector{Float64}(n_visited_states) # preallocate marginal vector
    @inbounds for i = 1:n_visited_states
        m[i] .= length(group_repeated_inds[i]) / npts
    end

    return m
end
