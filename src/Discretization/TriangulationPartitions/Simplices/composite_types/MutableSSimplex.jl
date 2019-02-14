import StaticArrays: MArray, MVector, SVector

###################################################################################
# Mutable simplex where vertices are represented by StaticArrays.MVector instances.
###################################################################################

"""
    MutableSSimplex(pts::Vector{SVector})

Construct a mutable `SSimplex` from a vector of `SVector`s.
"""
struct MutableSSimplex{dim, T} <: AbstractMutableSimplex{dim, T}
    vertices::Vector{MVector{dim, T}}

    function MutableSSimplex(pts::Vector{MVector{dim, T}}) where {dim, T}
        if !(length(pts) == dim + 1)
            err = """ The input cannot be converted to a simplex.
            Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
            """
            throw(DomainError(pts, err))
        end

        new{dim, T}(pts) 
    end
end

function Base.setindex!(s::MutableSSimplex, v::SVector{D, T}, i::Int) where {D, T}
    s[i] = v
end

function Base.setindex!(s::MutableSSimplex, v::SVector{D, T}, i::Int, j) where {D, T}
    s[i][j] = v
end


export MutableSSimplex
