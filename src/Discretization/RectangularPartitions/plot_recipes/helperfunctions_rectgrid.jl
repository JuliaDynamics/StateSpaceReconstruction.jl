import StaticArrays: SVector

#########################################################
# Helper functions to plot 3D rectangles.
#########################################################
"""
    rectangle3dvertices(x, y, z, ϵx, ϵy, ϵz)

Compute the coordinates of the vertices of a 3D rectangle, given the
coordinates of the origin (`x`, `y`, and `z`) and the corresponding edge lengths
(`ϵx`, `ϵy`, `ϵz`).
"""
function rectangle3dvertices(x, y, z, ϵx, ϵy, ϵz)
    v1b = [x,    y,    z]
    v2b = [x+ϵx, y,    z]
    v3b = [x+ϵx, y+ϵy, z]
    v4b = [x,    y+ϵy, z]
    v1t = [x,    y,    z+ϵz]
    v2t = [x+ϵx, y,    z+ϵz]
    v3t = [x+ϵx, y+ϵy, z+ϵz]
    v4t = [x,    y+ϵy, z+ϵz]
    (v1b, v2b, v3b, v4b, v1t, v2t, v3t, v4t)
end


"""
    rectangle3dvertices(o, ϵ)

Given the origin `o` (3-element vector) and edge lengths `ϵ` (also a
3-element vector) of a 3D rectangle, return the coordinates of its
vertices.
"""
function rectangle3dvertices(o, ϵ)
    x, y, z = (o...,)
    ϵx, ϵy, ϵz = (ϵ...,)
    rectangle3dvertices(x, y, z, ϵx, ϵy, ϵz)
end

"""
    connectvertices3d(o, ϵ)

Given the origin `o` (3-element vector) and edge lengths `ϵ` (also a
3-element vector) of a 3D rectangle, return a vector of line segments
(each a dim-by-n_vertices array) that when plotted together yields
the entire rectangle.
"""
function connectvertices3d(o, ϵ)
    # Get the vertices
    v1b, v2b, v3b, v4b, v1t, v2t, v3t, v4t = rectangle3dvertices(o, ϵ)

    # Connect vertices in top and bottom planes
    connect_top = zeros(3, 5)
    connect_bottom = zeros(3, 5)

    # Connect corners of the top and bottom planes
    corner1 = zeros(3, 2)
    corner2 = zeros(3, 2)
    corner3 = zeros(3, 2)
    corner4 = zeros(3, 2)
    connect_bottom[:, 1:5] = [v1b v2b v3b v4b v1b]
    connect_top[:, 1:5]    = [v1t v2t v3t v4t v1t]
    corner1[:, 1:2] = [v1b v1t]
    corner2[:, 1:2] = [v2b v2t]
    corner3[:, 1:2] = [v3b v3t]
    corner4[:, 1:2] = [v4b v4t]

    return [connect_top, connect_bottom,
            corner1, corner2, corner3, corner4]
end

export connectvertices3d


"""
    splitaxes(x)

Return a vector of the individual components of an array of points
provided as an array where each point is a column.
"""
splitaxes(x) = ([x[k, :] for k = 1:size(x, 1)]...,)

"""
    plot_3D_rect!(p, origin, edgelengths;
             lc = :black, lw = 1.0, ls = :solid)

Append a a 3D rectangle to an existing plot `p`, given the
origin `o` (3-element vector) and edge lengths
`ϵ` (3-element vector) of the rectangle.
"""
function plot_3D_rect!(p, o, ϵ;
             lc = :black, lw = 1.0, ls = :solid, lα = 0.7)

    linesegments = connectvertices3d(o, ϵ)
    for segment in linesegments
        plot!(p, splitaxes(segment),
            lc = lc, lw = lw, ls = ls, lα = lα)
    end
end

"""
Given the `origin` of a hyperrectangle and the `stepsizes` along each
coordinate axis, divide each axis into `N` stepsizes starting
from `origin[i] + (stepsizes[i]/N)/2`
"""
function fill_hyperrectangle_3D(origin, stepsizes, N = 10)
    xs, ys, zs = ([LinRange(origin[i] + stepsizes[i]/N/2, origin[i] + stepsizes[i] - stepsizes[i]/N/2, N) for i = 1:length(origin)]...,)
    pts = zeros(Float64, 3, length(xs)^3)
    i = 0
    for x in xs
        for y in ys
            for z in zs
                i += 1
                pts[:, i] = [x, y, z]
            end
        end
    end
    pts
end

struct Rectangle{D, T}
    origin::SVector{D, T}
    edgelengths::SVector{D, T}
end


rectangle3dvertices(r::Rectangle) = rectangle3dvertices((r.origin...,), (r.edgelengths...,))

function linesegments(r::Rectangle)
    # Get the vertices
    v1b, v2b, v3b, v4b, v1t, v2t, v3t, v4t = rectangle3dvertices(r.origin, r.edgelengths)

    # Connect vertices in top and bottom planes
    connect_top = zeros(3, 5)
    connect_bottom = zeros(3, 5)

    # Connect corners of the top and bottom planes
    corner1 = zeros(3, 2)
    corner2 = zeros(3, 2)
    corner3 = zeros(3, 2)
    corner4 = zeros(3, 2)
    connect_bottom[:, 1:5] = [v1b v2b v3b v4b v1b]
    connect_top[:, 1:5]    = [v1t v2t v3t v4t v1t]
    corner1[:, 1:2] = [v1b v1t]
    corner2[:, 1:2] = [v2b v2t]
    corner3[:, 1:2] = [v3b v3t]
    corner4[:, 1:2] = [v4b v4t]
    
    linesegments = [connect_top, connect_bottom, corner1, corner2, corner3, corner4]
    return linesegments
end

function connectrect(r::Rectangle) 
    v1b, v2b, v3b, v4b, v1t, v2t, v3t, v4t = rectangle3dvertices(r.origin, r.edgelengths)
    
    x = hcat(v1b, v2b, v2t, v3t, v3b, v4b, v4t, v1t, v1b, v1t, v2t, v2b, v3b, v3t, v4t, v4b, v1b)
    splitaxes(x)
end



export
Rectangle, 
connectrect, 
rectangle3dvertices, 
linesegments, 
connectvertices3d
