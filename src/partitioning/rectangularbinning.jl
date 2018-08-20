abstract type AbstractRectangularBinning <: Partition end

function Base.summary(eb::T) where T<:AbstractRectangularBinning
    npts = eb.n_pts
    dim = eb.dim
    binningtype = typeof(eb)
    return "$npts-point $dim-dimensional $binningtype"
end

function matstring(rb::T) where T<:AbstractRectangularBinning
    fields = fieldnames(rb)
    fields_str = String.(fields)
    maxlength = maximum([length(str) for str in fields_str]) + 2
    fields_str = [fields_str[i] *
                repeat(" ", maxlength - length(fields_str[i])) for i = 1:length(fields_str)]

    summaries = [join(":"*String(fields_str[i])*summary(getfield(rb, fields[i]))*"\n") for i = 1:length(fields_str)] |> join
    infoline = "The following fields are available:\n"

    return summary(rb)#*"\n\n"*infoline*summaries
end

Base.show(io::IO, rb::T) where {T<:AbstractRectangularBinning} = println(io, matstring(rb))

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
struct RectangularBinning <: AbstractRectangularBinning
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
struct EquidistantBinning <: AbstractRectangularBinning
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



dimension(rb::T) where T<:AbstractRectangularBinning = rb.dim
npoints(rb::T) where T<:AbstractRectangularBinning = rb.npts
lowerbound(rb::T) where T<:AbstractRectangularBinning = rb.bottom
upperbound(rb::T) where T<:AbstractRectangularBinning = rb.top
stepsizes(rb::T) where T<:AbstractRectangularBinning = rb.stepsizes
indices_nonempty_bins(rb::T) where T<:AbstractRectangularBinning = rb.inds_nonempty_bins
unique_nonempty_bins(rb::T) where T<:AbstractRectangularBinning = rb.unique_nonempty_bins
first_inds(rb::T) where T<:AbstractRectangularBinning = rb.first_inds
group_inds(rb::T) where T<:AbstractRectangularBinning = rb.group_inds
all_inds(rb::T) where T<:AbstractRectangularBinning = rb.all_inds

"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
function coords_bottom(b::AbstractRectangularBinning, i::Int)
    b.bottom .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end

"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
function coords_bottom(b::AbstractRectangularBinning, i::Int)
    b.bottom .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end


"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
coords_origin(b::AbstractRectangularBinning, i::Int) = coords_bottom(b, i)

"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
coords_origin(b::AbstractRectangularBinning, i::Int) = coords_bottom(b, i)


"""
Find the coordinates of the top of the i-th nonempty bin.
"""
function coords_top(b::AbstractRectangularBinning, i::Int)
    b.top .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end

"""
Find the coordinates of the top of the i-th nonempty bin.
"""
function coords_top(b::AbstractRectangularBinning, i::Int)
    b.top .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end

function unique_rows_info(embedding::Array{Float64, 2})
    first_inds = firstinds(groupslices(embedding, 1))
    group_inds = groupinds(groupslices(embedding, 1))
    all_inds = indexin_rows(embedding, unique(embedding, 1))
    first_inds, group_inds, all_inds
end


function unique_rows_info(A::Array{Int, 2})
    first_inds = firstinds(groupslices(A, 1))
    group_inds = groupinds(groupslices(A, 1))
    all_inds = indexin_rows(A, unique(A, 1))
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
    `bin_rectangular(E::T where T<:Embedding, stepsizes::Vector{Float64}) -> RectangularBinning`

Bin an `embedding` into rectangular boxes, specifying the bin size along each dimension
with `stepsizes`.
"""
function bin_rectangular(E::T, stepsizes::Vector{Float64}) where T<:AbstractEmbedding
    emb = transpose(E.points)
    n_pts, dim = size(emb, 2), size(emb, 1)

    bottom = [minimum(emb[:, i]) for i in 1:dim]
    top = [maximum(emb[:, i]) for i in 1:dim]
    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100

    # Indices of the bins. The coordinates of each point of the original
    # embedding are assigned an integer number indicating which bin along
    # the respective dimension it falls into.
    inds_nonempty_bins = zeros(Int, dim, n_pts)

    @inbounds for i = 1:n_pts
        inds_nonempty_bins[:, i] = ceil.(Int, (emb[:, i] .- bottom) ./ stepsizes)
    end

    slices = groupslices(inds_nonempty_bins, 2)
    first_inds = firstinds(slices)
    group_inds = groupinds(slices)

    unique_inds_nonempty_bins = unique(inds_nonempty_bins, 2)
    all_inds = radical_mad_improvement(inds_nonempty_bins, unique_inds_nonempty_bins)

    RectangularBinning(
                dim,
                n_pts,
                bottom,
                top,
                stepsizes,
                transpose(inds_nonempty_bins),
                transpose(unique_inds_nonempty_bins),
                first_inds,
                group_inds,
                all_inds)
end


function radical_mad_improvement(A, U)
    a = [A[:,i] for i in 1:size(A, 2)]
    u = [U[:,i] for i in 1:size(U, 2)]
    inds = Int[]
    npts = length(a)
    npts_unique = length(u)
    uvec = Vector{Int}(3)
    avec = Vector{Int}(3)

    nfound = 0
    for j = 1:npts
        for i = 1:npts_unique
            if a[j] == u[i]
                nfound += 1
                #print("""\tUnique point #$i found in pos $j.
                #        Setting the $nfound-th entry of inds to $i""")
                push!(inds, i)
            end
        end
    end
    inds
end

"""
    bin_rectangular(E::Embedding, n_bins::Int) -> RectangularBinning

Bin an embedding into `n_bins` rectangular, regularly sized boxes.
"""
function bin_rectangular(E::T, n_bins::Int) where T<:AbstractEmbedding
    emb = transpose(E.points)
    n_pts, dim = size(emb, 2), size(emb, 1)

    bottom = [minimum(emb[i, :]) for i in 1:dim]
    top = [maximum(emb[i, :]) for i in 1:dim]
    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100
    stepsizes = (top - bottom) / n_bins

    # Indices of the bins. The coordinates of each point of the original
    # embedding are assigned an integer number indicating which bin along
    # the respective dimension it falls into.
    inds_nonempty_bins = zeros(Int, dim, n_pts)

    @inbounds for i = 1:n_pts
        inds_nonempty_bins[:, i] = ceil.(Int, (emb[:, i] .- bottom) ./ stepsizes)
    end

    slices = groupslices(inds_nonempty_bins, 2)
    first_inds = firstinds(slices)
    group_inds = groupinds(slices)

    unique_inds_nonempty_bins = unique(inds_nonempty_bins, 2)
    all_inds = radical_mad_improvement(inds_nonempty_bins, unique_inds_nonempty_bins)

    bininfo = RectangularBinning(
                dim,
                n_pts,
                bottom,
                top,
                stepsizes,
                transpose(inds_nonempty_bins),
                transpose(unique_inds_nonempty_bins),
                first_inds,
                group_inds,
                all_inds)

    return bininfo
end


"""
    bin_rectangular(E::Embedding, boxsize_frac::Float64) -> RectangularBinning

Bin an `embedding` into rectangular boxes with regular lengths, with box dimensions given
as a fraction (0 < `boxsize_frac` < 1) of the data values along the associated coordinate
axis.
"""
function bin_rectangular(E::T, boxsize_frac::Float64) where T<:AbstractEmbedding
    emb = transpose(E.points)
    n_pts, dim = size(emb, 2), size(emb, 1)

    bottom = [minimum(emb[:, i]) for i in 1:dim]
    top = [maximum(emb[:, i]) for i in 1:dim]
    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100
    stepsizes = (top - bottom) * boxsize_frac

    # Indices of the bins. The coordinates of each point of the original
    # embedding are assigned an integer number indicating which bin along
    # the respective dimension it falls into.
    inds_nonempty_bins = zeros(Int, dim, n_pts)

    @inbounds for i = 1:n_pts
        inds_nonempty_bins[:, i] = ceil.(Int, (emb[:, i] .- bottom) ./ stepsizes)
    end

    slices = groupslices(inds_nonempty_bins, 2)
    first_inds = firstinds(slices)
    group_inds = groupinds(slices)

    unique_inds_nonempty_bins = unique(inds_nonempty_bins, 2)
    all_inds = radical_mad_improvement(inds_nonempty_bins, unique_inds_nonempty_bins)

    RectangularBinning(
                dim,
                n_pts,
                bottom,
                top,
                stepsizes,
                transpose(inds_nonempty_bins),
                transpose(unique_inds_nonempty_bins),
                first_inds,
                group_inds,
                all_inds)
end
