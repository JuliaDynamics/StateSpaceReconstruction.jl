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

unique_nonempty_bins::Array{Int, 2}
    Unique rows of inds_nonempty_bins.

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
struct RectangularBinning <: Partition
    dim::Int
    n_pts::Int
    bottom::Vector{Float64}
    top::Vector{Float64}
    stepsizes::Vector{Float64}
    inds_nonempty_bins::Array{Int, 2}
    unique_nonempty_bins::Array{Int, 2}
    first_inds::Vector{Int}
    group_inds::Vector{Vector{Int}}
    all_inds::Vector{Int}
end

""" An rectangular binning that has boxes of equal edge lengths. """
struct EquidistantBinning <:Partition
    dim::Int
    n_pts::Int
    bottom::Vector{Float64}
    top::Vector{Float64}
    stepsizes::Vector{Float64}
    inds_nonempty_bins::Array{Int, 2}
    unique_nonempty_bins::Array{Int, 2}
    first_inds::Vector{Int}
    group_inds::Vector{Vector{Int}}
    all_inds::Vector{Int}
end


dimension(rb::T where T<: RectangularBinning) = rb.dim
npoints(rb::T where T<: RectangularBinning) = rb.npts
lowerbound(rb::T where T<: RectangularBinning) = rb.bottom
upperbound(rb::T where T<: RectangularBinning) = rb.top
stepsizes(rb::T where T<: RectangularBinning) = rb.stepsizes
indices_nonempty_bins(rb::T where T<: RectangularBinning) = rb.inds_nonempty_bins
unique_nonempty_bins(rb::T where T<: RectangularBinning) = rb.unique_nonempty_bins
first_inds(rb::T where T<: RectangularBinning) = rb.first_inds
group_inds(rb::T where T<: RectangularBinning) = rb.group_inds
all_inds(rb::T where T<: RectangularBinning) = rb.all_inds

# TODO: rename all functions below to bin_rectangular (also in dependent libraries)


"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
function coords_bottom(b::EquidistantBinning, i::Int)
    b.bottom .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end

"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
function coords_bottom(b::RectangularBinning, i::Int)
    b.bottom .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end


"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
coords_origin(b::EquidistantBinning, i::Int) = coords_bottom(b, i)

"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
coords_origin(b::RectangularBinning, i::Int) = coords_bottom(b, i)


"""
Find the coordinates of the top of the i-th nonempty bin.
"""
function coords_top(b::EquidistantBinning, i::Int)
    b.top .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end

"""
Find the coordinates of the top of the i-th nonempty bin.
"""
function coords_top(b::RectangularBinning, i::Int)
    b.top .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end

function unique_rows_info(embedding::Array{Float64, 2})
    first_inds = firstinds(groupslices(embedding, 1))
    group_inds = groupinds(groupslices(embedding, 1))
    all_inds = indexin_rows(embedding, unique(embedding, 1))
    first_inds, group_inds, all_inds
end


function unique_rows_info(embedding::Array{Int, 2})
    first_inds = firstinds(groupslices(embedding, 1))
    group_inds = groupinds(groupslices(embedding, 1))
    all_inds = indexin_rows(embedding, unique(embedding, 1))
    first_inds, group_inds, all_inds
end

function indexin_rows(A1::Array{T, 2}, A2::Array{T, 2}) where {T<:Number}
    inds = Int[]
    for j = 1:size(A1, 1)
        for i = 1:size(A2, 1)
            if all(A1[j, :] .== A2[i, :])
                push!(inds, i)
            end
        end
    end
    return inds
end

"""
    `bin_rectangular(E::Embedding, boxsize_frac::Float64) -> RectangularBinning`

Bin an `embedding` into rectangular boxes with regular lengths, with box dimensions given
as a fraction (0< `boxsize_frac` < 1) of the data values along the associated coordinate
axis.
"""
function bin_rectangular(E::Embedding, boxsize_frac::Float64)
    emb = E.points
    n_pts, dim = size(emb, 1), size(emb, 2)

    bottom = [minimum(emb[:, i]) for i in 1:dim]
    top = [maximum(emb[:, i]) for i in 1:dim]
    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100
    stepsizes = (top - bottom) * boxsize_frac

    # Indices of the bins. The coordinates of each point of the original
    # embedding are assigned an integer number indicating which bin along
    # the respective dimension it falls into.
    inds_nonempty_bins = zeros{Int}(n_pts, dim)

    @inbounds for i = 1:n_pts
        inds_nonempty_bins[i, :] = ceil.(Int, (emb[i, :] - bottom) ./ stepsizes)
    end

    first_inds, group_inds, all_inds = unique_rows_info(inds_nonempty_bins)
    bininfo = RectangularBinning(
                dim,
                n_pts,
                bottom,
                top,
                stepsizes,
                inds_nonempty_bins,
                unique(inds_nonempty_bins, 1),
                first_inds,
                group_inds,
                all_inds)
end


"""
    bin_rectangular(E::Embedding, n_bins::Int) -> RectangularBinning

Bin an embedding into `n_bins` rectangular, regularly sized boxes.
"""
function bin_rectangular(E::T where T<:Embedding, n_bins::Int)
    emb = E.points
    n_pts, dim = size(emb, 1), size(emb, 2)

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
        inds_nonempty_bins[i, :] = ceil.(Int, (emb[i, :] - bottom) ./ stepsizes)
    end

    first_inds, group_inds, all_inds = unique_rows_info(inds_nonempty_bins)
    bininfo = RectangularBinning(
                dim,
                n_pts,
                bottom,
                top,
                stepsizes,
                inds_nonempty_bins,
                unique(inds_nonempty_bins, 1),
                first_inds,
                group_inds,
                all_inds)
end

"""
    bin_rectangular(E::Embedding, boxsize_frac::Float64) -> RectangularBinning

Bin an `embedding` into rectangular boxes with regular lengths, with box dimensions given
as a fraction (0 < `boxsize_frac` < 1) of the data values along the associated coordinate
axis.
"""
function bin_rectangular(E::T where T<:Embedding, boxsize_frac::Float64)
    emb = E.points
    n_pts, dim = size(emb, 1), size(emb, 2)

    bottom = [minimum(emb[:, i]) for i in 1:dim]
    top = [maximum(emb[:, i]) for i in 1:dim]
    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100
    stepsizes = (top - bottom) * boxsize_frac

    # Indices of the bins. The coordinates of each point of the original
    # embedding are assigned an integer number indicating which bin along
    # the respective dimension it falls into.
    inds_nonempty_bins = zeros{Int}(n_pts, dim)

    @inbounds for i = 1:n_pts
        inds_nonempty_bins[i, :] = ceil.(Int, (emb[i, :] - bottom) ./ stepsizes)
    end

    first_inds, group_inds, all_inds = unique_rows_info(inds_nonempty_bins)
    bininfo = RectangularBinning(
                dim,
                n_pts,
                bottom,
                top,
                stepsizes,
                inds_nonempty_bins,
                unique(inds_nonempty_bins, 1),
                first_inds,
                group_inds,
                all_inds)
end

"""
    `bin_rectangular(E::T where T<:Embedding, stepsizes::Vector{Float64}) -> RectangularBinning`

Bin an `embedding` into rectangular boxes, specifying the bin size along each dimension
with `stepsizes`.
"""
function bin_rectangular(E::T where T<:Embedding, stepsizes::Vector{Float64})
    emb = E.points
    n_pts, dim = size(emb, 1), size(emb, 2)

    bottom = [minimum(emb[:, i]) for i in 1:dim]
    top = [maximum(emb[:, i]) for i in 1:dim]
    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100

    # Indices of the bins. The coordinates of each point of the original
    # embedding are assigned an integer number indicating which bin along
    # the respective dimension it falls into.
    inds_nonempty_bins = zeros(Int, n_pts, dim)

    @inbounds for i = 1:n_pts
        inds_nonempty_bins[i, :] = ceil.(Int, (emb[i, :] - bottom) ./ stepsizes)
    end

    first_inds, group_inds, all_inds = unique_rows_info(inds_nonempty_bins)
    bininfo = RectangularBinning(
                dim,
                n_pts,
                bottom,
                top,
                stepsizes,
                inds_nonempty_bins,
                unique(inds_nonempty_bins, 1),
                first_inds,
                group_inds,
                all_inds)
end
