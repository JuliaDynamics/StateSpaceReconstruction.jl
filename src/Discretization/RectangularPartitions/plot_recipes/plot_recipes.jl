#include("helperfunctions_rectgrid.jl")

using RecipesBase

@recipe function f(rp::RectangularPartition)
    label --> ""
    for rect in rp.rectangles
        @series begin
            connectrect(rect)
        end
    end
end

@recipe function f(r::Rectangle)
    label --> ""
    lc --> :black
    connectrect(r)
end

