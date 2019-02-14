using RecipesBase

@recipe function f(s::T) where {T <: AbstractSimplex}
    seriestype := :path
    splitaxes(connectvertices(s))
end


@recipe function f(simplices::Vararg{T,N}) where {T <: AbstractSimplex, N}
    seriestype := :path
    legend --> false

    for simplex in simplices
        @series begin
            label --> ""
            splitaxes(connectvertices(simplex))
        end
    end

end


@recipe function f(simplices::AbstractVector{T}) where {T <: AbstractSimplex}
    seriestype := :path
    legend --> false

    for simplex in simplices
        @series begin
            label --> ""
            splitaxes(connectvertices(simplex))
        end
    end

end
