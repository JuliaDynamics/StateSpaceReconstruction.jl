
export
point_contained,
point_contained!


"""
point_contained(simplex::ST, point::Vector{VT}) where {T <: AbstractSimplex, VT}

Check if `point` is contained in the interior of `simplex`.

The `point` is contained within `simplex` if it can be expressed as a
convex combination of the vertices of `simplex`. To check this condition,
we iteratively check the signed volume of a modified simplex. Specifically,
at the i-th iteration step, substitute the i-th vertex of `simplex` with
`point` and take the signed volume when when one vertex is substituted.
In the next step, replace the substitated vertex by the original, then
replace the (i + 1)th vertex by `point` and repeat the process.
If any two consecutive signs are different, the point is not contained.

When the point is very close to the edge of a simplex, numerical accuracies
will play a role, so the function sometimes (very rarely) puts `point` inside
`simplex`.
"""
function point_contained(simplex::T, point::Vector{VT}) where {T <: AbstractSimplex, VT}
    dim = dimension(simplex)

    # The signs of the linear combination coefficients
    coefficient_signs = zeros(Float64, dim + 1)

    # A temporary simplex into which we substitute `point`
    temp_simplex = MutableSimplex([zeros(Float64, dim) for i = 1:(dim + 1)])

    # Using the temporary simplex, replace the the first vertex by `point`
    # and fill the remaining vertices with the original `simplex` vertices.
    temp_simplex[1] = point
    [temp_simplex[i] = simplex[i] for i = 2:(dim + 1)]

    # The signed volume of the simplex
    coefficient_signs[1] = sign(orientation(temp_simplex))

    # Iteratively replace the other vertices and see if the signed
    # volume changes. It it does, the point can not be expressed as a
    # convex combination of the vertices of `simplex`.
    for κ = 2:dim

        # Replace the i-th vertex of the temporary simplex with `point`
        temp_simplex[κ] = point

        # Replace the remaining vertices of the temporary simplex
        # with the original vertices
        remaining_vertex_idxs = setdiff(1:(dim + 1), κ)
        [temp_simplex[i] = simplex[i] for i in remaining_vertex_idxs]

        # Check the signed volume and stop if it has changed.
        coefficient_signs[κ] = sign(orientation(temp_simplex))

        if !(coefficient_signs[κ-1] == coefficient_signs[κ])
            return false
        end
    end

    # Check when replacing the last point.
    temp_simplex[end] = point
    [temp_simplex[i] = simplex[i] for i = 1:dim]

    coefficient_signs[end] = sign(orientation(temp_simplex))

    if !(coefficient_signs[end-1] == coefficient_signs[end])
        return false
    else
        return true
    end
end

"""
    point_contained!(temp_simplex::MutableSimplex, simplex::ST, point::Vector{VT}) where {T <: AbstractSimplex, VT}

Check if `point` is contained in the interior of `simplex`.

The `point` is contained within `simplex` if it can be expressed as a
convex combination of the vertices of `simplex`. To check this condition,
we iteratively check the signed volume of a modified simplex. Specifically,
at the i-th iteration step, substitute the i-th vertex of `simplex` with
`point` and take the signed volume when when one vertex is substituted.
In the next step, replace the substitated vertex by the original, then
replace the (i + 1)th vertex by `point` and repeat the process.
If any two consecutive signs are different, the point is not contained.

When the point is very close to the edge of a simplex, numerical accuracies
will play a role, so the function sometimes (very rarely) puts `point` inside
`simplex`.
"""
function point_contained!(temp_simplex::MutableSimplex, simplex::T, point::Vector{VT}) where {T <: AbstractSimplex, VT}
    dim = dimension(simplex)

    # The signs of the linear combination coefficients
    coefficient_signs = zeros(Float64, dim + 1)

    # Using the temporary simplex, replace the the first vertex by `point`
    # and fill the remaining vertices with the original `simplex` vertices.
    temp_simplex[1] = point
    [temp_simplex[i] = simplex[i] for i = 2:(dim + 1)]

    # The signed volume of the simplex
    coefficient_signs[1] = sign(orientation(temp_simplex))

    # Iteratively replace the other vertices and see if the signed
    # volume changes. It it does, the point can not be expressed as a
    # convex combination of the vertices of `simplex`.
    for κ = 2:dim

        # Replace the i-th vertex of the temporary simplex with `point`
        temp_simplex[κ] = point

        # Replace the remaining vertices of the temporary simplex
        # with the original vertices
        remaining_vertex_idxs = setdiff(1:(dim + 1), κ)
        [temp_simplex[i] = simplex[i] for i in remaining_vertex_idxs]

        # Check the signed volume and stop if it has changed.
        coefficient_signs[κ] = sign(orientation(temp_simplex))

        if !(coefficient_signs[κ-1] == coefficient_signs[κ])
            return false
        end
    end

    # Check when replacing the last point.
    temp_simplex[end] = point
    [temp_simplex[i] = simplex[i] for i = 1:dim]

    coefficient_signs[end] = sign(orientation(temp_simplex))

    if !(coefficient_signs[end-1] == coefficient_signs[end])
        return false
    else
        return true
    end
end
