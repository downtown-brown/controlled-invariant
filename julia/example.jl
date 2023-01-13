using Polyhedra
using IntervalArithmetic

function example()
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

    S = I∞(Ω, U, ϵ, A, B, phi, psi)
end

function inverted_pendulum(case=1, ϵ=0.001)
    n = 2
    dt = 0.01

    m = 0.2
    g = 9.8
    l = 0.3
    J = 0.006
    b = 0.1

    # TODO
    A(x) = [1 dt; m*g*l/J*dt*cos(x[1]) (1 - dt*b/J)]
    B(x)::Matrix{Float64} = reshape(dt*[0; l/J*cos(x[1])], 2, 1)

    phi(X, x) = IntervalBox(dt*[0; sin(X[1]) - cos(x[1])*x[1] + sin(x[1])])
    psi(X, U, x) = dt*U[1]*IntervalBox([0; l/J*(cos(X[1]) - cos(x[1]))])

    if case == 1
        X = IntervalBox(-0.05..0.05, -0.01..0.01)
        U = (-0.1..0.1)
    elseif case == 2
        X = IntervalBox(0.10..0.17, -0.01..0.01)
        U = (-0.4..(-0.1))
    elseif case == 3
        ϵ = 10.
        X = IntervalBox(-170..170, -180..180)
        U = (-1..1)
    end

    U_p = i2p(IntervalBox(U))

    Ω = [X]

    S = I∞(Ω, U, ϵ, A, B, phi, psi)

    plot(legend=:none)
    for s in S
        plot!(s)
    end

    return S
end


function cpre(Ω::Vector{IntervalBox{M,T}},
              U, ϵ, A, B, phi, psi) where {M,T<:Real}

    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}()
    E = Vector{IntervalBox{M,T}}()
    L = Vector{IntervalBox{M,T}}(Ω)

    U_p = i2p(IntervalBox(U))

    Ω_p = [i2p(o) for o in merge_best_effort(Ω)]

    i = 0;
    i_s = i_n = i_e = 0
    while length(L) != 0
        x = popfirst!(L)
        #println("considering interval ", x)

        x_p = i2p(x)

        x_m = mid(x)
        image_poly = A(x_m)*x_p + i2p(phi(x, x_m)) + i2p(psi(x, U, x_m))
        image_over = image_poly + B(x_m)*U_p

        Xi = image_over .∩ Ω_p

        if !intersects(image_over, Ω_p)
            #println("rejected")
            push!(N, x)
        elseif can_translate_into(image_poly, Xi)
            #println("accepted")
            push!(S, x)
        elseif diam(x) < ϵ
            #println("boundary")
            push!(E, x)
        else
            #println("bisected")
            x₁, x₂ = bisect(x, 0.5)
            push!(L, x₁)
            push!(L, x₂)
        end
        i += 1;
    end

    println("Considered ", i, " intervals")

    return (S, N, E)
end


function I∞(Ω::Vector{IntervalBox{M,T}},
            U, ϵ::T, A, B, phi, psi) where {M,T<:Real}

    Y = Vector{IntervalBox{M,T}}()
    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}([emptyinterval() × emptyinterval()])
    E = Vector{IntervalBox{M,T}}([emptyinterval() × emptyinterval()])
    i = 0
    while !isempty(N) || !isempty(E)
        #Y = merge_best_effort(Ω)
        Y = Ω
        (S, N, E) = cpre(Y, U, ϵ, A, B, phi, psi)
        Ω = S
        i += 1
    end

    return merge_best_effort(S)

end

function i2p(A::IntervalBox{M,T}) where {M,T<:Real}
    n = length(A)
    res = Vector{HalfSpace{Float64, Vector{Float64}}}(undef, 2*n)
    for i in 1:n
        tmp = zeros(n)
        tmp[i] = -1
        res[i] = HalfSpace(tmp, -A.v[i].lo)
        tmp = zeros(n)
        tmp[i] = 1
        res[n+i] = HalfSpace(tmp, A.v[i].hi)
    end

    return polyhedron(hrep(res))
end

function intersects(A, B)

    for B_i in B
        if volume(A ∩ B_i) > 1e-16
            return true
        end
    end

    return false
end
