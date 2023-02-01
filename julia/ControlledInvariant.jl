import Base.merge, Base.intersect
using IntervalArithmetic
import IntervalArithmetic.⊂
using Polyhedra

using Plots

using CDDLib
lib = CDDLib.Library(:float)

#lib = DefaultLibrary{Float64}()

################################################################################
## Example functions with continuous input
################################################################################
function example_synthetic()
    n = 2;
    mu = 0.5;
    dt = 0.01;


    X = IntervalBox(-1.0..1.0,-1.0..1.0);
    U = (-2.0..2.0);
    U_p = i2p(IntervalBox(U))

    A(x) = [1 dt; dt 1];
    phi(X, x) = 0;
    B(x)::Matrix{Float64} = reshape(dt*[(1-mu)*x[1]; 4*(1-mu)*x[2]], 2, 1);
    psi(X,U,x) = U[1]*dt*IntervalBox([(1-mu)*(X[1]-x[1]); 4*(1-mu)*(X[2]-x[2])]);

    ϵ = .0001;

    Ω = [X]

    S = I_infty(Ω, U, ϵ, A, B, phi, psi)
end

p = nothing
function example_pendulum_direct(case=1, ϵ=0.001)
    global p
    n = 2
    dt = 0.01

    m = 0.2
    g = 9.8
    l = 0.3
    J = 0.006
    b = 0.1

    A(x) = [1 dt; m*g*l/J*dt*cos(x[1]) (1 - dt*b/J)]
    B(x)::Matrix{Float64} = reshape(dt*[0; l/J*cos(x[1])], 2, 1)

    phi(X, x) = IntervalBox([0; m*g*l/J*dt*(sin(X[1]) - cos(x[1])*x[1])])
    psi(X, U, x) = dt*U[1]*IntervalBox([0; l/J*(cos(X[1]) - cos(x[1]))])

    if case == 1
        X = IntervalBox(-0.05..0.05, -0.01..0.01)
        U = IntervalBox(-0.1..0.1)
    elseif case == 2
        X = IntervalBox(-(pi)..(pi), -1.0..1.0)
        U = IntervalBox(-0.1..0.1)
    elseif case == 3
        X = IntervalBox(-180..180, -180..180)
        U = IntervalBox(-1..1)
    end

    p = ProblemData(A, phi, B, psi, U, ϵ)

    Ω = Vector{IntervalData}([IntervalData(X)])

    S, N, E = I_direct(Ω, p)

    return S
end


function example_pendulum(case=1, ϵ=0.001)
    n = 2
    dt = 0.01

    m = 0.2
    g = 9.8
    l = 0.3
    J = 0.006
    b = 0.1

    A(x) = [1 dt; m*g*l/J*dt*cos(x[1]) (1 - dt*b/J)]
    B(x)::Matrix{Float64} = reshape(dt*[0; l/J*cos(x[1])], 2, 1)

    phi(X, x) = IntervalBox([0; m*g*l/J*dt*(sin(X[1]) - cos(x[1])*x[1])])
    psi(X, U, x) = dt*U[1]*IntervalBox([0; l/J*(cos(X[1]) - cos(x[1]))])

    if case == 1
        X = IntervalBox(-0.05..0.05, -0.01..0.01)
        U = (-0.1..0.1)
    elseif case == 2
        X = IntervalBox(-(pi)..(pi), -1.0..1.0)
        U = (-0.1..0.1)
    elseif case == 3
        X = IntervalBox(-180..180, -180..180)
        U = (-1..1)
    end

    Ω = [X]

    S = I_infty(Ω, U, ϵ, A, B, phi, psi)

    return S
end

function inverted_pendulum(ϵ=0.001)
    n = 2
    dt = 0.001

    m = 1
    g = 9.81
    l = 1
    J = m*l^2

    A(x) = [1 dt; m*g*l/J*dt*cos(x[1]) 1]
    B(x)::Matrix{Float64} = reshape(dt*[0; 1/J], 2, 1)

    phi(X, x) = IntervalBox(dt*[0.; m*g*l/J*(sin(X[1]) - cos(x[1])*x[1])])
    psi(X, U, x) = IntervalBox([0.; 0.])

    X = IntervalBox(-180.0..180.0, -180.0..180.0)
    U = (-5.0..5.0)

    U_p = i2p(IntervalBox(U))

    Ω = [X]

    S, Us = I_infty(Ω, U, ϵ, A, B, phi, psi)

    return S, Us
end


################################################################################
## Example functions for switched systems
################################################################################
function example_pendulum_sw(ϵ=0.001)
    n = 2
    dt = 0.01

    m = 0.2
    g = 9.8
    l = 0.3
    J = 0.006
    b = 0.1

    X = IntervalBox(-0.05..0.05, -0.01..0.01)

    X = IntervalBox(-(pi)..(pi), -(pi)..(pi))

    U = range(-0.5,stop=0.5,step=0.05)

    f = Vector{Function}()
    for u in U
        tmp(x) = IntervalBox([
            x[1] + dt*x[2];
            (1-dt*b/J)*x[2] + dt*(m*g*l/J*sin(x[1]) + l/J*cos(x[1])*u)
        ])
        push!(f, tmp)
    end

    Ω = [X]
    Ω_c = [emptyinterval() × emptyinterval()]

    I_infty_sw(Ω, Ω_c, f, ϵ)
end

function example_dcdc()
    x_c = 70
    r_c = 0.005
    x_l = 3
    r_l = 0.05
    r_0 = 1
    v_s = 1

    A_1 = [-r_l/x_l 0;
           0 -1/(x_c*(r_c+r_0))]
    A_2 = [-1/x_l*(r_l + r_0*r_c/(r_c + r_0)) -r_0/(x_l*(r_c+r_0));
           r_0/(x_c*(r_c+r_0)) -1/(x_c*(r_c+r_0))]
    b = [v_s/x_l; 0]

    t_s = 0.5

    A_1 = exp(A_1*t_s)
    A_2 = exp(A_2*t_s)

    b_1 = [0.165974147222481;
           0]

    b_2 = [0.165872918752656;
           0.000589016663477]

    f_1(x) = A_1*x + b_1
    f_2(x) = A_2*x + b_2
    f = [f_1, f_2]
    Ω = [(1.15..1.55) × (1.09..1.17)]
    Ω_c = [emptyinterval() × emptyinterval()]
    ϵ = 0.001

    I_infty_sw(Ω, Ω_c, f, ϵ)
end

function inverted_pendulum_sw(ϵ=0.001)
    n = 2
    dt = 0.01

    m = 1
    g = 9.81
    l = 1
    J = m*l^2

    U = (-5.0..5.0)
    U = range(-5.0,stop=5.0,step=1.)

    f = Vector{Function}()
    for u in U
        tmp(x) = IntervalBox([
            x[1] + dt*x[2];
            x[2] + dt*(m*g*l/J*sin(x[1]) + u/J)
        ])
        push!(f, tmp)
    end

    X = IntervalBox(-(pi)..(pi), -(pi)..(pi))

    U_p = i2p(IntervalBox(U))

    Ω = [X]
    Ω_c = [emptyinterval() × emptyinterval()]

    S = I_infty_sw(Ω, Ω_c, f, ϵ,)

    return S
end

################################################################################
## Controlled invariant functions with continuous input
################################################################################
struct ProblemData2
    A::Function
    Phi::Function
    B::Function
    Psi::Function
    U::IntervalBox
    U_p::Polyhedron
    ϵ::Float64
    ProblemData2(
        A::Function,
        Phi::Function,
        B::Function,
        Psi::Function,
        U::IntervalBox,
        ϵ::Float64
    ) = new(A, Phi, B, Psi, U, i2p(U), ϵ)
end
ProblemData = ProblemData2

mutable struct IntervalData5
    interval::IntervalBox
    poly::Polyhedron
    lchild::Ref{IntervalData5}
    rchild::Ref{IntervalData5}
    P_u_over::Polyhedron
    P_over::Polyhedron
    invariant::Bool
    iter::Int64
    function IntervalData5(x::IntervalBox, p::ProblemData)
        poly = i2p(x)
        x_m = mid(x)
        P_u_over = p.A(x_m)*poly + i2p(p.Phi(x, x_m)) + i2p(p.Psi(x, p.U, x_m))
        P_over = P_u_over + p.B(x_m)*p.U_p
        new(x, poly, Ref{IntervalData5}(), Ref{IntervalData5}(), P_u_over, P_over, true, 0)
    end

end
IntervalData = IntervalData5

xdbg = nothing
sdbg = []
ndbg = []
edbg = []
function I_direct(Ω::Vector{IntervalData}, p::ProblemData)
    global xdbg, sdbg, edbg, ndbg
    S = Vector{IntervalData}()
    N = Vector{IntervalData}()
    E = Vector{IntervalData}()
    L = Vector{IntervalData}(Ω)


    i = 0;

    k = 1;
    Ω_p = [i2p(o) for o in merge_best_effort([n.interval for n in L])]
    Nc = reduce(convexhull, Ω_p)
    Nd = regiondiff(Nc, Ω_p)


    while length(L) != 0
        x = popfirst!(L)
        if x.iter == k
            break
        end
        xdbg = x
        sdbg = S
        ndbg = N
        edbg = E

        if !intersects(x.P_over, Ω_p)
            push!(N, x)
            x.invariant = false
            k += 1
            S = Vector{IntervalData}()

            Ω_p = [i2p(o) for o in merge_best_effort([n.interval for n in L])]
            Nc = reduce(convexhull, Ω_p)
            Nd = regiondiff(Nc, Ω_p)

        elseif can_translate_into(x.P_u_over, x.P_over, Nc, Nd)
            push!(S, x)
            x.iter = k
            push!(L, x)

        elseif diam(x.interval) < p.ϵ
            push!(E, x)
            x.invariant = false
            k += 1
            S = Vector{IntervalData}()

            Ω_p = [i2p(o) for o in merge_best_effort([n.interval for n in L])]
            Nc = reduce(convexhull, Ω_p)
            Nd = regiondiff(Nc, Ω_p)

        else
            x₁, x₂ = bisect(x, p)
            prepend!(L, [x₁])
            prepend!(L, [x₂])
        end
        i += 1;
        if mod(i, 100) == 0
            println("considered ", i, " intervals")
        end
    end

    println("considered ", i, " intervals")


    return (S, N, E)
end

function IntervalArithmetic.bisect(x::IntervalData, p::ProblemData)
    l, r = bisect(x.interval, 0.5)
    x_l = IntervalData(l, p)
    x.lchild = x_l
    x_r = IntervalData(r, p)
    x.rchild = x_r
    return (x_l, x_r)
end

k = 1
"""
Compute an approximation of the one step backward reachable set of Ω, intersected with Ω
"""
a = []
b = []
c = []
cs = []
ds = []
es = []
cr = []
dr = []
er = []
function I_approx(Ω::Vector{IntervalBox{M,T}},
              U, ϵ, A, B, phi, psi) where {M,T<:Real}
    global k, a, b, c

    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}()
    E = Vector{IntervalBox{M,T}}()
    L = Vector{IntervalBox{M,T}}(Ω)

    U_p = i2p(IntervalBox(U))

    Ω_p = [i2p(o) for o in merge_best_effort(Ω)]

    i = 0;
    i_s = i_n = i_e = 0

    Nc = reduce(convexhull, Ω_p)
    Nd = regiondiff(Nc, Ω_p)

    push!(a, Ω_p)
    push!(b, Nc)
    push!(c, Nd)

    while length(L) != 0
        x = popfirst!(L)

        x_p = i2p(x)

        x_m = mid(x)
        image_poly = A(x_m)*x_p + i2p(phi(x, x_m)) + i2p(psi(x, U, x_m))
        image_over = image_poly + B(x_m)*U_p

        if !intersects(image_over, Ω_p)
            push!(N, x)

            if k == 2
                push!(cs, x)
                push!(ds,image_poly)
                push!(es,image_over)
            end

        elseif can_translate_into(image_poly, image_over, Nc, Nd)
            push!(S, x)
            if k == 2
                push!(cr, x)
                push!(dr,image_poly)
                push!(er,image_over)
            end

        elseif diam(x) < ϵ
            push!(E, x)
        else
            x₁, x₂ = bisect(x, 0.5)
            push!(L, x₁)
            push!(L, x₂)
        end
        i += 1;
    end

    println("considered ", i, " intervals")

    f = open("test.txt", "a")
    print_interval_list(f, S, "S", k)
    print_interval_list(f, N, "N", k)
    print_interval_list(f, E, "E", k)
    println(f, "")
    close(f)

    k += 1

    return (S, N, E)
end

function I_final(Ω::Vector{IntervalBox{M,T}},
              U, ϵ, A, B, phi, psi) where {M,T<:Real}
    global k

    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}()
    E = Vector{IntervalBox{M,T}}()
    L = Vector{IntervalBox{M,T}}(Ω)

    U_p = i2p(IntervalBox(U))

    Ω_p = [i2p(o) for o in merge_best_effort(Ω)]

    i = 0;
    i_s = i_n = i_e = 0

    Nc = reduce(convexhull, Ω_p)
    Nd = regiondiff(Nc, Ω_p)

    Ux = []

    while length(L) != 0
        x = popfirst!(L)

        x_p = i2p(x)

        x_m = mid(x)
        image_poly = A(x_m)*x_p + i2p(phi(x, x_m)) + i2p(psi(x, U, x_m))
        image_over = image_poly + B(x_m)*U_p

        Xi = image_over .∩ Ω_p

        Us = get_translate_into(image_poly, image_over, Nc, Nd)
        for u in Us
            v = collect(points(vrep(u)))
            if length(v) != 2
                println("Error")
                continue
            end
            if v[1][2] >= v[2][2]
                U_x = IntervalBox([x[1].lo, x[2].lo, v[2][2]/B(x_m)[2]],
                                  [x[1].hi, x[2].hi, v[1][2]/B(x_m)[2]])
            else
                U_x = IntervalBox([x[1].lo, x[2].lo, v[1][2]/B(x_m)[2]],
                                  [x[1].hi, x[2].hi, v[2][2]/B(x_m)[2]])

            end
            push!(Ux, U_x)
        end
    end

    return Ux
end

function print_interval_list(f, S, name, i)
    println(f, name, "{", i, "} = {")
    for s in S
        println(f, "interval([", s[1].lo, ";", s[2].lo, "],[", s[1].hi, ";", s[2].hi, "]),")
    end
    println(f, "}")
end

function I_approx_par(Ω::Vector{IntervalBox{M,T}},
              U, ϵ, A, B, phi, psi) where {M,T<:Real}

    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}()
    E = Vector{IntervalBox{M,T}}()
    L = Vector{IntervalBox{M,T}}(Ω)

    U_p = i2p(IntervalBox(U))

    Ω_p = [i2p(o) for o in Ω]

    i = 0;
    i_s = i_n = i_e = 0

    Nc = reduce(convexhull, Ω_p)
    Nd = regiondiff(Nc, Ω_p)

    while length(L) != 0
        x = popfirst!(L)

        x_p = i2p(x)

        x_m = mid(x)
        image_poly = A(x_m)*x_p + i2p(phi(x, x_m)) + i2p(psi(x, U, x_m))
        image_over = image_poly + B(x_m)*U_p

        Xi = image_over .∩ Ω_p

        if !intersects(image_over, Ω_p)
            push!(N, x)
        elseif can_translate_into(image_poly, image_over, Nc, Nd)
            push!(S, x)
        elseif diam(x) < ϵ
            push!(E, x)
        else
            x₁, x₂ = bisect(x, 0.5)
            push!(L, x₁)
            push!(L, x₂)
        end
        i += 1;
    end

    println("Considered ", i, " intervals")

    return (S, N, E)
end

"""
Compute an approximation of the limit of the I operator
"""
function I_infty(Ω::Vector{IntervalBox{M,T}},
            U, ϵ::T, A, B, phi, psi) where {M,T<:Real}

    f = open("test.txt", "w")
    close(f)
    Y = Vector{IntervalBox{M,T}}()
    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}([emptyinterval() × emptyinterval()])
    E = Vector{IntervalBox{M,T}}([emptyinterval() × emptyinterval()])
    i = 0
    while (!isempty(N) || !isempty(E)) && !isempty(Ω)
        #Y = merge_best_effort(Ω)
        Y = Ω
        (S, N, E) = I_approx(Y, U, ϵ, A, B, phi, psi)
        Ω = S
        i += 1
    end

    #Us = I_final(Ω, U, ϵ, A, B, phi, psi)

    return merge_best_effort(S)#, Us

end

################################################################################
## Controlled invariant functions with continuous input
################################################################################

"""
Compute an approximation of the I operator for switched systems
"""
kk = 1
function I_sw(Ω::Vector{IntervalBox{M,T}},
              Ω_c::Vector{IntervalBox{M,T}},
              X::IntervalBox{M,T},
              f::Vector{Function},
              ϵ::T) where {M,T<:Real}

    global kk

    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}()
    E = Vector{IntervalBox{M,T}}()
    L = Vector{IntervalBox{M,T}}(Ω)

    i = 0;
    i_s = i_n = i_e = 0
    while length(L) != 0
        x = pop!(L)

        images = [f_p(x) for f_p in f]
        if any([(!intersects(image, Ω_c) && (image ⊂ X )) for image in images])
            push!(S, x)
        elseif all([!intersects(image, Ω) for image in images])
            push!(N, x)
        elseif diam(x) < ϵ
            push!(E, x)
        else
            x₁, x₂ = bisect(x, 0.5)
            push!(L, x₁)
            push!(L, x₂)
        end
        i += 1;
    end

    println("Considered ", i, " intervals")


    f = open("sw.txt", "a")
    print_interval_list(f, S, "S", kk)
    print_interval_list(f, N, "N", kk)
    print_interval_list(f, E, "E", kk)
    println(f, "")
    close(f)

    kk += 1

    return (S, N, E)
end

"""
Compute an approximation of the limit of the I operator
"""
function I_infty_sw(Ω::Vector{IntervalBox{M,T}},
            Ω_c::Vector{IntervalBox{M,T}},
            f::Vector{Function},
            ϵ::T) where {M,T<:Real}

    X = IntervalBox(Ω[1])
    Y_c = Vector{IntervalBox{M,T}}(Ω_c)
    Y = Vector{IntervalBox{M,T}}()
    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}([emptyinterval() × emptyinterval()])
    E = Vector{IntervalBox{M,T}}([emptyinterval() × emptyinterval()])
    i = 0


    f1 = open("sw.txt", "w")
    close(f1)
    while !isempty(N) || !isempty(E)
        Y = Ω
        Y_c = merge_best_effort([Y_c; N; E])
        (S, N, E) = I_sw(Y, Y_c, X, f, ϵ)
        Ω = S
        i += 1
    end

    return (S, Y_c)

end


################################################################################
## Interval and polyhedral functions
################################################################################

"""
Convert an interval to polyhedron
"""
function i2p(A::IntervalBox{M,T}) where {M,T<:Real}
    n = length(A)
    res = Vector{HalfSpace{T, Vector{T}}}(undef, 2*n)
    for i in 1:n
        tmp = zeros(T, n)
        tmp[i] = -1
        res[i] = HalfSpace(tmp, -A.v[i].lo)
        tmp = zeros(T, n)
        tmp[i] = 1
        res[n+i] = HalfSpace(tmp, A.v[i].hi)
    end

    return polyhedron(hrep(res), lib)
end

"""
Check if `A` ∩ `B` is nonempty for polyhedron `A` and union of polyhedra `B`
"""
function intersects(A::Polyhedron, B::Vector{<:Polyhedron})
    for B_i in B
        if dim(A ∩ B_i) != 0
            return true
        end
    end

    return false
end

"""
Check if `A` ∩ `B` is empty for interval `A` and union of intervals `B`
"""
function intersects(A::IntervalBox{M,T},
                    B::Vector{IntervalBox{M,T}}) where {M,T<:Real}

    for B_i in B
        if !isempty(A ∩ B_i)
            return true
        end
    end

    return false
end


function test_polyhedra()
    C = polyhedron(vrep([[4, 4],[2, 4], [1, 5]]), lib);

    P_v = [
        polyhedron(vrep([[2., 2.], [1., 1.], [-5., 3.]]), lib);
        polyhedron(vrep([[2., 2.], [1., 1.], [5., -3.]]), lib);
        polyhedron(vrep([[-3., -4.], [1., 1.], [-5., 3.]]), lib);
        polyhedron(vrep([[2., 2.], [6., 3.], [5., -3.]]), lib);
    ]

    P = polyhedron(vrep(([[-5., 3.], [-3., -4.], [5., -3.], [6.,3.]])), lib);
end

"""
Compute the set { u | C + u ⊂ N } for polyhedra `C` and `N`
"""
function translate_into(C::Polyhedron, N::Polyhedron)
    res = Vector{HalfSpace{Float64, Vector{Float64}}}()
    for n in halfspaces(hrep(N))
        max = -Inf
        for c in points(vrep(C))
            tmp = n.a'*c
            if tmp > max
                max = tmp
            end
        end
        push!(res, HalfSpace(n.a, n.β - max))
    end

    return polyhedron(hrep(res), lib)
end

"""
Compute the set { u | C + u ⊂ X \\ N } for polyhedra `C` and `N` and union of polyhedra `N`
"""
function can_translate_into(C::Polyhedron, C_over::Polyhedron, Nc::Polyhedron, Nd::Vector{<:Polyhedron})

    Nc = C_over ∩ Nc
    Ndd = Vector{Polyhedron}()
    for n in Nd
        tmp = C_over ∩ n
        if dim(tmp) >= 0
            push!(Ndd, tmp)
        end
    end

    U1 = translate_into(C, Nc)
    U2 = translate_touching(C, Ndd)

    return !(U1 ⊂ U2)
end

function get_translate_into(C::Polyhedron, C_over::Polyhedron, Nc::Polyhedron, Nd::Vector{<:Polyhedron})

    Nc = C_over ∩ Nc
    Ndd = Vector{Polyhedron}()
    for n in Nd
        tmp = C_over .∩ n
        if dim(tmp) >= 0
            push!(Ndd, tmp)
        end
    end


    U1 = translate_into(C, Nc)
    U2 = translate_touching(C, Ndd)

    return regiondiff(U1, U2)
end

"""
Compute the set { u | C + u ∩ N ≠ ∅ } for polyhedra `C` and `N`
"""
function translate_touching(C::Polyhedron, N::Polyhedron)::Polyhedron

    res = Vector{HalfSpace{Float64, Vector{Float64}}}()

    for n in halfspaces(hrep(N))
        min = Inf
        for c in points(vrep(C))
            tmp = n.a'*c
            if tmp < min
                min = tmp
            end
        end
        push!(res, HalfSpace(n.a, n.β - min))
    end

    try
        for c in halfspaces(hrep(C))
            min = Inf
            for n in points(vrep(N))
                tmp = c.a'*n
                if tmp < min
                    min = tmp
                end
            end
            push!(res, HalfSpace(-c.a, c.β - min))
        end
    catch
        println(C)
        error()
    end

    return polyhedron(hrep(res), lib)
end

"""
Compute the set { u | C + u ≠ ∅ } for polyhedron `C` and union of polyhedra `N`
"""
function translate_touching(C::Polyhedron, N::Vector{<:Polyhedron})::Vector{Polyhedron}
    return [translate_touching(C, n) for n in N]
end

Pdbg = Nothing
function ⊂(P::Polyhedron, Q::Vector{T})::Bool where {T<:Polyhedron}
    global Pdbg
    Pdbg = P
    if dim(P) < 0
        return true
    end
    if isempty(Q)
        return false
    end

    P = copy(P)

    nQ = length(Q)
    k = 1

    while dim(P ∩ Q[k]) < 0
        k += 1
        if k > nQ
            return false
        end
    end

    for q_h in halfspaces(hrep(Q[k]))
        S = copy(P)

        d1 = dim(S)

        intersect!(S, HalfSpace(-q_h.a, -q_h.β))

        if dim(S) < d1
            continue
        end

        if k == nQ
            return false
        elseif !(S ⊂ Q[k+1:end])
            return false
        end

        intersect!(P, q_h)
    end

    return true
end

"""
Computes the set difference `P` \\ `Q`between a
single polyhedron `P` and a union of polyhedra `Q`
"""
function regiondiff(P::T, Q::Vector{T})where {T<:Polyhedron}
    res = Vector{T}()

    if isempty(Q)
        return res
    end

    P = copy(P)

    nQ = length(Q)

    k = 1
    while dim(P ∩ Q[k]) == 0
        k += 1
        if k > nQ
            return [P]
        end
    end

    for q_h in halfspaces(hrep(Q[k]))
        S = copy(P)

        intersect!(S, HalfSpace(-q_h.a, -q_h.β))

        if dim(S) < 2
            continue
        end

        if k < nQ
            append!(res, regiondiff(S, Q[k+1:end]))
        else
            push!(res, S)
        end

        intersect!(P, q_h)
    end

    return res
end



function merge_best_effort(A::Vector{T}) where T<:IntervalBox
    results = A
    n = length(A)
    while true
        results = merge_all_dims(results)
        n == length(results) && break
        n = length(results)
    end
    results
end

function merge_all_dims(A::Vector{T}) where T<:IntervalBox
    results = A
    if length(A) > 1
        for d in 1:length(A[1])
            cmp(x, y) = lt_d(x, y, d)
            results = merge(results, cmp)
        end
    end
    results
end


function merge(A::Vector{T}, B::Vector{T}, comp=lt) where T<:IntervalBox
    merged = Vector{T}()
    n = length(A) + length(B)
    i = 1
    j = 1
    while i <= length(A) || j <= length(B)
        if i <= length(A) && j<= length(B)
            combo = merge(A[i], B[j])
            if isa(combo, T)
                victim = combo
                i = i + 1
                j = j + 1
            elseif comp(A[i], B[j])
                victim = A[i]
                i = i + 1
            else
                victim = B[j]
                j = j + 1
            end
        elseif i <= length(A)
            victim = A[i]
            i = i + 1
        else
            victim = B[j]
            j = j + 1
        end

        if length(merged) > 0
            combo = merge(victim, merged[end])
            if isa(combo, T)
                merged[end] = combo
            else
                push!(merged, victim)
            end
        else
            push!(merged, victim)
        end
    end
    return merged
end

function merge(A::IntervalBox{M, T}, B::IntervalBox{M, T}) where {M, T<:Real}
    if A ⊆ B || B ⊆ A
        return A ∪ B
    end
    mustbeequal = false
    for d in 1:M
        if A[d] ∩ B[d] == ∅
            return (A, B)
        elseif A[d] != B[d]
            if mustbeequal
                return (A, B)
            else
                mustbeequal = true
            end
        end
    end

    return A ∪ B
end

function merge(A::Vector{T}, comp=lt) where T<:IntervalBox

    if length(A) > 1

        halfway = div(length(A), 2)
        left = A[1:halfway]
        right = A[halfway+1:end]

        left = merge(left)
        right = merge(right)

        return merge(left, right)
    else
        return A
    end
end

function lt(A::IntervalBox{M, T}, B::IntervalBox{M, T}) where {M, T<:Real}
    for d in 1:M
        if A[d] != B[d] && A[d].lo > B[d].lo
            return false
        end
    end
    return true
end
