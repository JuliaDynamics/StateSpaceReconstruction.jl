
import Simplices.Delaunay.delaunay

"""
    DelaunayTriangulation

A Delaunay triangulation in dimension D. If `d`
is an instance of `DelaunayTriangulation`, then
`d.indices[i]` gives the D + 1 indices of the vertices
corresponding to the i-th simplex. The indices are
expressed in terms of the points it was produced
from.
"""
struct DelaunayTriangulation
    indices::AbstractArray{Int32, 2}
end


function delaunaytriang(d::Dataset)
    triang = delaunay(transpose(Matrix(d)))
    DelaunayTriangulation(hcat(triang...))
end


export delaunaytriang, DelaunayTriangulation
