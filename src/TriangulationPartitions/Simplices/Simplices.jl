using Reexport

@reexport module Simplices

include("AbstractSimplex.jl")

struct Simplex{T} <: AbstractSimplex
    vertices::Vector{Vector{T}}
end


function Simplex(pts)
    s = size(pts)
    if (maximum(s) - minimum(s)) > 1
        throw(DomainError(pts, "The input cannot be converted to a simplex. size(pts) must be (dim, dim + 1) or (dim + 1, dim)"))
    end

    if s[1] > s[2] # Rows are points
        return Simplex([pts[i, :] for i = 1:maximum(s)])
    end

    if s[2] > s[1] # Columns are points
        return Simplex([pts[:, i] for i = 1:maximum(s)])
    end

end

export Simplex

end # module
