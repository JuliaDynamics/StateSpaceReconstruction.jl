"""
Equidistant rectangular binning of an embedding.

dim::Int
    Dimension of the space

n_pts::Int
    How many points does the binning cover?

bottom::Vector{Float64}
    Coordinates of the origin/bottom corner of the binning

top::Vector{Float64}
    Coordinates of the top corner of the binning.

stepsizes::Vector{Float64}
    Edge lengths for the bins along each axis.

inds_nonempty_bins::Array{Int, 2}
    Indices of nonempty bins. One bin per row.

first_inds::Vector{Int}
    Row indices of `inds_nonempty_bins` for which the unique rows of `inds_nonempty_bins`
    first appears.

group_inds::Vector{Vector{Int}}
    For each unique row in `inds_nonempty_bins`, a vector of row indices indicating where
    in `inds_nonempty_bins` rows equal to that unique row is located.

all_inds::Vector{Int}

    `all_inds` has the same number of rows as `inds_nonempty_bins`. The ith index of
    `all_inds` indicates which row of `unique(inds_nonempty_bins)` corresponds to that
    particular bin.
"""

@with_kw struct EquidistantBinning <: Partition
    dim::Int = 0
    n_pts::Int = 0
    bottom::Vector{Float64} = Float64[]
    top::Vector{Float64} = Float64[]
    stepsizes::Vector{Float64} = Float64[]
    inds_nonempty_bins::Array{Int, 2} = Array{Int, 2}(0, 0)
    first_inds::Vector{Int} = Int[]
    group_inds::Vector{Vector{Int}} = Vector{Vector{Int}}(0)
    all_inds::Vector{Int} = Int[]
end


"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
function coords_bottom(b::EquidistantBinning, i::Int)
    b.bottom .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end

"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
coords_origin(b::EquidistantBinning, i::Int) = coords_bottom(b, i)

"""
Find the coordinates of the top of the i-th nonempty bin.
"""
function coords_top(b::EquidistantBinning, i::Int)
    b.top .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end


function unique_rows_info(embedding)
    first_inds = firstinds(groupslices(embedding, 1))
    group_inds = groupinds(groupslices(embedding, 1))
    all_inds = indexin_rows(embedding, unique(embedding, 1))
    first_inds, group_inds, Int.(all_inds)
end

function indexin_rows(A1::Array{Float64, 2}, A2::Array{Float64, 2})
    inds = []
    for j = 1:size(A1, 1)
        for i = 1:size(A2, 1)
            if all(A1[j, :] .== A2[i, :])
                push!(inds, i)
            end
        end
    end
    return inds
end

function indexin_rows(A1::Array{Int, 2}, A2::Array{Int, 2})
    inds = []
    for j = 1:size(A1, 1)
        for i = 1:size(A2, 1)
            if all(A1[j, :] .== A2[i, :])
                push!(inds, i)
            end
        end
    end
    return inds
end

function bin_equidistant(embedding::Embedding, n_bins::Int)
    emb = embedding.points
    dim = embedding.dim
    n_pts = size(emb, 1)

    bottom = [minimum(emb[:, i]) for i in 1:dim]
    top = [maximum(emb[:, i]) for i in 1:dim]
    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100
    stepsizes = (top - bottom) / n_bins

    # Indices of the bins. The coordinates of each point of the original
    # embedding are assigned an integer number indicating which bin along
    # the respective dimension it falls into.
    inds_nonempty_bins = zeros(Int, n_pts, dim)

    @inbounds for i = 1:n_pts
        inds_nonempty_bins[i, :] = ceil((emb[i, :] - bottom) ./ stepsizes)
    end

    first_inds, group_inds, all_inds = unique_rows_info(inds_nonempty_bins)
    bininfo = EquidistantBinning(
                dim = dim,
                n_pts = n_pts,
                bottom = bottom, top = top, stepsizes = stepsizes,
                inds_nonempty_bins = inds_nonempty_bins,
                first_inds = first_inds,
                group_inds = group_inds,
                all_inds = all_inds)
end

bin_equidistant(pts::AbstractArray{Float64, 2}, n_bins::Int) = bin_equidistant(Embedding(pts), b_bins)
bin_equidistant(pts::AbstractArray{Int, 2}, n_bins::Int) = bin_equidistant(float.(pts), b_bins)
