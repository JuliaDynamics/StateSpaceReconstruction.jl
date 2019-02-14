using RecipesBase

@recipe function f(triang::AbstractTriangulationPartition) where {D, T}
    simplices =
    legend --> false
    gridalpha --> 0.5
    foreground_color_grid --> :black
    for point in getpoints(triang)
        @series begin
            seriestype := :scatter
            ms --> 2,
            mc --> :black
            fc --> :black

            splitaxes(point)
        end
    end
    for simplex in getsimplices(triang)
        @series begin
            lc --> :black
            lw --> 0.8
            lα --> 0.8
            simplex
        end
    end
end
