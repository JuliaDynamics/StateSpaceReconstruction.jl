import StaticArrays: SArray, SVector

##################################################################################
# Static simplex where vertices are represented by StaticArrays.SVector instances.
##################################################################################

struct SSimplex{dim, T} <: AbstractImmutableSimplex{dim, T}
    vertices::Vector{SVector{dim, T}}

    function SSimplex(pts::Vector{SVector{dim, T}}) where {dim, T}
        if !(length(pts) == dim + 1)
            err = """ The input cannot be converted to a simplex.
            Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
            """
            throw(DomainError(pts, err))
        end

        new{dim, T}(pts) 
    end

    function SSimplex(pts::Vector{Vector{T}}) where {T}
        dim = length(pts[1])
        if !(length(pts) == dim + 1)
            err = """ The input cannot be converted to a simplex.
            Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
            """
            throw(DomainError(pts, err))
        end

        [SVector{dim, T}(pt) for pt in pts]
    end
end

export SSimplex
