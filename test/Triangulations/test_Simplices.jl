
using StateSpaceReconstruction
using Test

# Constructors
@test Simplex([rand(3) for i = 1:4]) isa Simplex
@test Simplex(rand(3, 4)) isa Simplex
@test Simplex(rand(4, 3)) isa Simplex
@test Simplex(rand(4, 5)) isa Simplex
@test Simplex(rand(5, 4)) isa Simplex
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
