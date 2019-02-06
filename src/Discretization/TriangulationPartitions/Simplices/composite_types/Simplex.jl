import StaticArrays:
	SArray,
    MArray

################################################################################
# Static simplex where vertices are represented by some type of abstract vector.
################################################################################

struct Simplex{D, T} <: AbstractImmutableSimplex{D, T}
    vertices::Vector{AbstractVector{T}}
end

"""
    Simplex(pts::Vector{AbstractVector})

Construct a simplex from a vector of vectors.
"""
function Simplex(pts::Vector{Vector{T}}) where {T}
    if !(length(pts) == length(pts[1]) + 1)
        err = """ The input cannot be converted to a simplex.
        Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
        """
        throw(DomainError(pts, err))
    end
    dim = length(pts[1])
    Simplex{dim, T}([pts[i] for i = 1:length(pts)])
end

"""
    Simplex(pts::Vector{SVector{D, T}})

Construct a simplex from a vector of static vectors.
"""
function Simplex(pts::Vector{SArray{Tuple{D},T,1,D}}) where {D, T}
    if !(length(pts) == length(pts[1]) + 1)
        err = """ The input cannot be converted to a simplex.
        Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
        """
        throw(DomainError(pts, err))
    end
    dim = length(pts[1])
    Simplex{dim, T}([pts[i] for i = 1:length(pts)])
end

"""
    Simplex(pts::Vector{MVector{D, T}})

Construct a simplex from a vector of mutable static vectors.
"""
function Simplex(pts::Vector{MArray{Tuple{D},T,1,D}}) where {D, T}
    if !(length(pts) == length(pts[1]) + 1)
        err = """ The input cannot be converted to a simplex.
        Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
        """
        throw(DomainError(pts, err))
    end
    dim = length(pts[1])
    Simplex{dim, T}([pts[i] for i = 1:length(pts)])
end

"""
    Simplex(pts::AbstractArray{T, 2}) where T

Construct a simplex from an array of size `(dim + 1)-by-(dim)` or
`(dim)-by-(dim + 1)` (faster).
"""
function Simplex(pts::AbstractArray{T, 2}) where {T}
    s = size(pts)

    if (maximum(s) - minimum(s)) > 1
        err = """ The input cannot be converted to a simplex.
            size(pts) must be (dim, dim + 1) or (dim + 1, dim).
        """
        throw(DomainError(pts, err))
    end

    # Assuming rows of the input are vertices
    if s[1] > s[2]
        dim = s[2]
        return Simplex{dim, T}([pts[i, :] for i = 1:maximum(s)])
    end

    # Assuming columns of the input are vertices
    if s[2] > s[1]
        dim = s[1]
        return Simplex{dim, T}([pts[:, i] for i = 1:maximum(s)])
    end
end



export Simplex
