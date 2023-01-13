using Polyhedra
using Plots

C = polyhedron(vrep([[4, 4],[2, 4], [1, 5]]));

P_v = [
    polyhedron(vrep([[2., 2.], [1., 1.], [-5., 3.]]));
    polyhedron(vrep([[2., 2.], [1., 1.], [5., -3.]]));
    polyhedron(vrep([[-3., -4.], [1., 1.], [-5., 3.]]));
    polyhedron(vrep([[2., 2.], [6., 3.], [5., -3.]]));
]

P = polyhedron(vrep(([[-5., 3.], [-3., -4.], [5., -3.], [6.,3.]])));

"""
Compute the set { u | C + u ⊂ N } for polyhedra `C` and `N`
"""
function translate_into(C::Polyhedron, N::Polyhedron)
    res = Vector{HalfSpace{Float64, Vector{Float64}}}()
    for n in hrep(N).halfspaces
        max = -Inf
        for c in vrep(C).points.points
            tmp = n.a'*c
            if tmp > max
                max = tmp
            end
        end
        push!(res, HalfSpace(n.a, n.β - max))
    end

    return polyhedron(hrep(res))
end

"""
Compute the set { u | C + u ⊂ N } for polyhedron `C` and union of polyhedra `N`
"""
function translate_into(C, N)

    Nc = P
    Nd = regiondiff(Nc, N)

    U1 = translate_into(C, Nc)
    U2 = translate_touching(C, Nd)

    return regiondiff(U1, U2)
end

"""
Compute the set { u | C + u ⊂ X \\ N } for polyhedra `C` and `N` and union of polyhedra `N`
"""
function can_translate_into(C::Polyhedron, N::Vector{<:Polyhedron})
    Nc = reduce(convexhull, N)
    Nd = regiondiff(Nc, N)

    U1 = translate_into(C, Nc)
    U2 = translate_touching(C, Nd)

    return !(U1 ⊂ U2)
end

"""
Compute the set { u | C + u ∩ N ≠ ∅ } for polyhedra `C` and `N`
"""
function translate_touching(C::Polyhedron, N::Polyhedron)

    res = Vector{HalfSpace{Float64, Vector{Float64}}}()

    for n in hrep(N).halfspaces
        min = Inf
        for c in vrep(C).points.points
            tmp = n.a'*c
            if tmp < min
                min = tmp
            end
        end
        push!(res, HalfSpace(n.a, n.β - min))
    end

    for c in hrep(C).halfspaces
        min = Inf
        for n in vrep(N).points.points
            tmp = c.a'*n
            if tmp < min
                min = tmp
            end
        end
        push!(res, HalfSpace(-c.a, c.β - min))
    end

    return polyhedron(hrep(res))
end

"""
Compute the set { u | C + u ≠ ∅ } for polyhedron `C` and union of polyhedra `N`
"""
function translate_touching(C::Polyhedron, N::Vector{<:Polyhedron})
    return [translate_touching(C, n) for n in N]
end

"""
Computes the set difference `P \\ Q` between two single polyhedra
"""
function Base.:\(P::Polyhedron, Q::Polyhedron)
    if isempty(P ∩ Q)
        return [P]
    end

    res = []

    for i in 1:length(Q.constraints)
        tmp = copy(P)
        addconstraint!(tmp, LazySets.HalfSpace(-Q.constraints[i].a,
                                               -Q.constraints[i].b))

        for j in 1:i-1
            addconstraint!(tmp, Q.constraints[j])
        end
        push!(res, tmp)

    end

    return res
end

"""
Computes the set difference `P` \\ `Q`between a
single polyhedron `P` and a union of polyhedra `Q`
"""
function Base.:\(P::Polyhedron, Q::Vector{T}) where {T<:Polyhedron}

    S = [P]

    for q in Q
        for s in S
            append!(res, s \ q)
        end
        S = copy(res)
    end

    return res
end

function subset_nonrecursive(P::Polyhedron, Q::Vector{T})::Bool where {T<:Polyhedron}
    res = Vector{T}()
    P = copy(P)


    nQ = length(Q)
    k = 1
    while volume(P ∩ Q[k]) < 1e-16
        k += 1
        if k > nQ
            return false
        end
    end

    q_h = hrep(Q[k]).halfspaces

    for i in 1:length(q_h)
        s_h = copy(hrep(P).halfspaces)
        push!(s_h, HalfSpace(-q_h[i].a, -q_h[i].β))

        S = polyhedron(hrep(s_h))

        if volume(S) < 1e-16
            continue
        end

        if k == nQ
            return false
        elseif !(S ⊂ Q[k+1:end])
            return false
        end

        push!(hrep(P).halfspaces, q_h[i])
    end

    return true
end

function ⊂(P::Polyhedron, Q::Vector{T})::Bool where {T<:Polyhedron}
    if isempty(vrep(P).points)
        return true
    end
    if isempty(Q)
        return false
    end

    P = copy(P)

    nQ = length(Q)
    k = 1
    while volume(P ∩ Q[k]) < 1e-16
        k += 1
        if k > nQ
            return false
        end
    end

    q_h = hrep(Q[k]).halfspaces

    for i in 1:length(q_h)
        s_h = copy(hrep(P).halfspaces)
        push!(s_h, HalfSpace(-q_h[i].a, -q_h[i].β))

        S = polyhedron(hrep(s_h))

        if volume(S) < 1e-16
            continue
        end

        if k == nQ
            return false
        elseif !(S ⊂ Q[k+1:end])
            return false
        end

        push!(hrep(P).halfspaces, q_h[i])
    end

    return true
end

i = 0
function regiondiff(P::T, Q::Vector{T})where {T<:Polyhedron}
    global i

    res = Vector{T}()

    if isempty(Q)
        return res
    end

    P = copy(P)

    nQ = length(Q)

    k = 1
    while volume(P ∩ Q[k]) < 1e-16
        k += 1
        if k > nQ
            return [P]
        end
    end

    q_h = hrep(Q[k]).halfspaces

    for i in 1:length(q_h)
        s_h = copy(hrep(P).halfspaces)
        push!(s_h, HalfSpace(-q_h[i].a, -q_h[i].β))

        S = polyhedron(hrep(s_h))

        if volume(S) < 1e-16
            continue
        end

        if k < nQ
            append!(res, regiondiff(S, Q[k+1:end]))
        else
            push!(res, S)
        end

        push!(hrep(P).halfspaces, q_h[i])
    end

    i += 1

    return res
end


struct regiondiffstate2
    P::Polyhedron
    Q::Vector{<:Polyhedron}
    parent::regiondiffstate2
    processed::Bool
end

function regiondiff_nonrecursive(P::Polyhedron, Q::Vector{T})where {T<:Polyhedron}
    global i
    res = Vector{T}()

    if isempty(Q)
        return res
    end

    P = copy(P)

    nQ = length(Q)
    k = 1
    while volume(P ∩ Q[k]) < 1e-16
        k += 1
        if k > nQ
            return [P]
        end
    end

    q_h = hrep(Q[k]).halfspaces

    for i in 1:length(q_h)
        s_h = copy(hrep(P).halfspaces)
        push!(s_h, HalfSpace(-q_h[i].a, -q_h[i].β))

        S = polyhedron(hrep(s_h))

        if volume(S) < 1e-16
            continue
        end

        if k < nQ
            append!(res, regiondiff(S, Q[k+1:end]))
        else
            push!(res, S)
        end

        push!(hrep(P).halfspaces, q_h[i])
    end

    return res
end
