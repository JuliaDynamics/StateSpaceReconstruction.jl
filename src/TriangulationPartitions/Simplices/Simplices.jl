using Reexport

@reexport module Simplices

include("AbstractSimplex.jl")
import StaticArrays.SArray

struct Simplex{T} <: AbstractSimplex{T}
    vertices::Vector{Vector{T}}
end

mutable struct MutableSimplex{T} <: AbstractSimplex{T}
    vertices::Vector{Vector{T}}
end

"""
    Simplex(pts::AbstractArray{T, 2}) where T

Construct a simplex from a matrix.
"""
function Simplex(pts::AbstractArray{T, 2}) where T
    s = size(pts)
    if (maximum(s) - minimum(s)) > 1
        throw(DomainError(pts, "The input cannot be converted to a simplex. size(pts) must be (dim, dim + 1) or (dim + 1, dim)"))
    end

    if s[1] > s[2] # Rows are points
        return Simplex{T}([pts[i, :] for i = 1:maximum(s)])
    end

    if s[2] > s[1] # Columns are points
        return Simplex{T}([pts[:, i] for i = 1:maximum(s)])
    end
end

"""
    Simplex(pts::Vector{AbstractVector})

Construct a simplex from a vector of vectors.
"""
function Simplex(pts::Vector{AbstractVector{T}}) where T
    if !(length(pts) == length(pts[1]) + 1)
        throw(DomainError(pts, "The input cannot be converted to a simplex. Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices."))
    end
    Simplex{T}([pts[i] for i = 1:length(pts)])
end

"""
    Simplex(pts::Vector{SVector})

Construct a simplex from a vector of `SVector`s.
"""
function Simplex(pts::Vector{SArray{Tuple{D}, T, 1, D}}) where {T, D}
    if !(length(pts) == D + 1)
        throw(DomainError(pts, "The input cannot be converted to a simplex. Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices."))
    end
    Simplex{T}([pts[i] for i = 1:length(pts)])
end

function MutableSimplex(pts::AbstractArray{T, 2}) where T
    s = size(pts)
    if (maximum(s) - minimum(s)) > 1
        throw(DomainError(pts, "The input cannot be converted to a simplex. size(pts) must be (dim, dim + 1) or (dim + 1, dim)"))
    end

    if s[1] > s[2] # Rows are points
        return MutableSimplex([pts[i, :] for i = 1:maximum(s)])
    end

    if s[2] > s[1] # Columns are points
        return MutableSimplex([pts[:, i] for i = 1:maximum(s)])
    end
end

"""
    Simplex(pts::Vector{AbstractVector})

Construct a simplex from a vector of vectors.
"""
function MutableSimplex(pts::Vector{AbstractVector{T}}) where T
    if !(length(pts) == length(pts[1]) + 1)
        throw(DomainError(pts, "The input cannot be converted to a simplex. Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices."))
    end
    MutableSimplex{T}([pts[i] for i = 1:length(pts)])
end

"""
    MutableSimplex(pts::Vector{SVector})

Construct a simplex from a vector of `SVector`s.
"""
function MutableSimplex(pts::Vector{SArray{Tuple{D}, T, 1, D}}) where {T, D}
    if !(length(pts) == D + 1)
        throw(DomainError(pts, "The input cannot be converted to a simplex. Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices."))
    end
    MutableSimplex{T}([pts[i] for i = 1:length(pts)])
end

export
Simplex,
MutableSimplex

end # module
