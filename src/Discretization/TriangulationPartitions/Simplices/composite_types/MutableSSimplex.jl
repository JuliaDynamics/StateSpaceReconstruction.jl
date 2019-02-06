import StaticArrays: MArray

###################################################################################
# Mutable simplex where vertices are represented by StaticArrays.MVector instances.
###################################################################################

mutable struct MutableSSimplex{D, T} <: AbstractMutableSimplex{D, T}
    vertices::Vector{MArray{Tuple{D}, T, 1, D}}
end


function Base.setindex!(s::MutableSSimplex, v::SArray{Tuple{D}, T, 1, D},
        i::Int) where {D, T}
    s[i] = v
end

function Base.setindex!(s::MutableSSimplex, v::SArray{Tuple{D}, T, 1, D},
        i::Int, j) where {D, T}
    s[i][j] = v
end


# """
#     MutableSSimplex(pts::Vector{SVector})

# Construct a mutable `SSimplex` from a vector of `SVector`s.
# """
# function MutableSSimplex(pts::Vector{SArray{Tuple{D}, T, 1, D}}) where {D, T}
#     if !(length(pts) == D + 1)
#         err = """ The input cannot be converted to a simplex.
#         Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
#         """
#         throw(DomainError(pts, err))
#     end
#     MutableSSimplex{D, T}([pts[i] for i = 1:length(pts)])
# end


export MutableSSimplex
