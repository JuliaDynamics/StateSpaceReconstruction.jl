
import Test
import StaticArrays:
    SVector,
    MVector
using StateSpaceReconstruction: Simplices


# Constructors
@testset "Immutable simplices" begin
    @test Simplex([rand(3) for i = 1:4]) isa Simplex
    @test Simplex(rand(3, 4)) isa Simplex
    @test Simplex(rand(4, 3)) isa Simplex
    @test Simplex(rand(4, 5)) isa Simplex
    @test Simplex(rand(5, 4)) isa Simplex
    @test SSimplex([SVector{3}(rand(3)) for i = 1:4]) isa SSimplex

    @test_throws DomainError Simplex(rand(5, 3))

    s = Simplex([rand(3) for i = 1:4])

    @test npoints(s) == 4
    @test dimension(s) == 3

    @test volume(s) isa Float64
    @test volume(s) >= 0

    @test orientation(s) isa Float64


    s1 = Simplex(rand(3, 4))
    s2 = Simplex(rand(4, 3))

    # Intersection volume
    @test intersect(s1, s2) >= 0

    # Intersecting vertices
    @test intersect_vertices(s1, s2) isa Array{Float64, 2}
end



# Constructors
@testset "Mutable simplices" begin
    @test MutableSimplex([rand(3) for i = 1:4]) isa MutableSimplex
    @test MutableSimplex([rand(3) for i = 1:4]) isa MutableSimplex
    @test MutableSSimplex([MVector{3}(rand(3)) for i = 1:4]) isa MutableSSimplex

    @test MutableSimplex(rand(3, 4)) isa MutableSimplex
    @test MutableSimplex(rand(4, 3)) isa MutableSimplex
    @test MutableSimplex(rand(4, 5)) isa MutableSimplex
    @test MutableSimplex(rand(5, 4)) isa MutableSimplex
    @test_throws DomainError MutableSimplex(rand(5, 3))

    s = MutableSimplex([rand(3) for i = 1:4])

    @test volume(s) isa Float64
    @test volume(s) >= 0

    @test orientation(s) isa Float64


    s1 = MutableSimplex(rand(3, 4))
    s2 = MutableSimplex(rand(4, 3))

    # Intersection volume
    @test intersect(s1, s2) >= 0

    # Intersecting vertices
    @test intersect_vertices(s1, s2) isa Array{Float64, 2}
end

@testset "Volume intersections" begin
    # Generate the vertices of two intersecting simplices, then construct simplices of
    # different types from those vertices. All types should give the same intersecting
    # volume and intersecting vertices.
    vertices_s1 = [rand(3) for i = 1:4]

    s1 = Simplex(vertices_s1)
    s1_mutable = MutableSimplex(vertices_s1)
    #ss1 = SSimplex([SVector{3, Float64}(vertices_s1[i]) for i = 1:4])
    ss1_mutable = MutableSSimplex([MVector{3, Float64}(vertices_s1[i]) for i = 1:4])


    s2 = generate_intersecting_simplex(s1)
    vertices_s2 = s2.vertices

    s2_mutable = MutableSimplex(vertices_s2)
    ss2 = SSimplex([SVector{3, Float64}(vertices_s2[i]) for i = 1:4])
    ss2_mutable = MutableSSimplex([MVector{3, Float64}(vertices_s2[i]) for i = 1:4])

    # Individual intersections should all be nonzero
    @test s1 ∩ s2 > 1e-12
    @test s1 ∩ s2_mutable > 1e-12
    @test s1 ∩ ss2 > 1e-12
    @test s1 ∩ ss2_mutable > 1e-12

    # Using the different types should give the same answer. Just try a few different
    # combinations.
    @test s1 ∩ s2 == s1 ∩ s2
    @test s1 ∩ s2 == s1 ∩ ss2
    @test s1 ∩ s2 == s1 ∩ s2_mutable
    @test s1 ∩ s2 == s1 ∩ ss2_mutable

end


@testset "Does a point lie inside simplex?" begin
    s = Simplex(rand(3, 4))
    n = 10000
    contained = Vector{Bool}(undef, n)

    # Interior points
    for i = 1:n
        p = generate_interior_points(s, 1)[1]
        contained[i] = point_contained(s, p)
    end

    @test all(x -> (x .== true), contained)

    # Exterior points
    for i = 1:n
        p = generate_exterior_points(s, 1)[1]
        contained[i] = point_contained(s, p)
    end

    @test all(x -> (x .== false), contained)
end


@testset "Refining simplices" begin
    
    pts = [rand(3) for i = 1:20]
    inv_pts = invariantize(pts);
    triang = triangulate(inv_pts[1:end - 1])
    triang_im = triangulate(inv_pts[2:end])

    simplices = [Simplex(inv_pts[triang.simplexindices[i]]) for i in 1:length(triang.simplexindices)];
    simplices_im = [Simplex(inv_pts[triang_im.simplexindices[i]]) for i in 1:length(triang_im.simplexindices)];
   
    @test refine(simplices[1]) isa Vector{Simplex{3,Float64}}
    @test refine(simplices[1], simplices[2]) isa Tuple{Vector{Simplex{3,Float64}}, Vector{Simplex{3,Float64}}}

    k = 2
    dim = 3
    rules = simplex_splitting_rules(k, dim)
    @test refine(simplices[1], k = 2, splitting_rules = rules) isa Vector{Simplex{3,Float64}}
    @test refine(simplices[1], simplices[2], k = 2, splitting_rules = rules) isa Tuple{Vector{Simplex{3,Float64}}, Vector{Simplex{3,Float64}}} 
end


@testset "Subsampling" begin
    dim = 3
    inv_pts = invariantize([rand(dim) for i = 1:20])
    triang = triangulate(inv_pts[1:end-1])
    triang_im = triangulate(inv_pts[2:end])
    simplices = [Simplex(inv_pts[triang.simplexindices[i]]) for i in 1:length(triang.simplexindices)];
    simplices_im = [Simplex(inv_pts[triang_im.simplexindices[i]]) for i in 1:length(triang_im.simplexindices)];
    
    coeffs = subsample_coeffs(dim, 50, false)

    # On single simplices.
    @test subsample(simplices[1], coeffs) isa Vector{SVector{3, Float64}}
    @test subsample(simplices[1]) isa Vector{SVector{3, Float64}}

    # Pairwise
    ss_simplices_coeffs_precomputed = subsample(simplices[1], simplices_im[1], coeffs)
    ss_simplices_coeffs_computed_onthefly = subsample(simplices[1], simplices_im[1])
    @test ss_simplices_coeffs_precomputed isa Tuple{Vector{SVector{3, Float64}}, Vector{SVector{3, Float64}}}
    @test ss_simplices_coeffs_computed_onthefly isa Tuple{Vector{SVector{3, Float64}}, Vector{SVector{3, Float64}}}
end