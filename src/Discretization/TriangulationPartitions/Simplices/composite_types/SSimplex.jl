import StaticArrays: SArray

##################################################################################
# Static simplex where vertices are represented by StaticArrays.SVector instances.
##################################################################################

struct SSimplex{D, T} <: AbstractImmutableSimplex{D, T}
    vertices::Vector{SArray{Tuple{D}, T, 1, D}}
end


# """
#     SSimplex(pts::Vector{SVector})

# Construct an `SSimplex` from a vector of `SVector`s.
# """
# function SSimplex(pts::Vector{SArray{Tuple{D}, T, 1, D}}) where {T, D}
#     if !(length(pts) == D + 1)
#         err = """ The input cannot be converted to a simplex.
#         Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
#         """
#         throw(DomainError(pts, err))
#     end

#     SSimplex{D, T}([pts[i] for i = 1:length(pts)])
# end


export SSimplex
