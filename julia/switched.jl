using IntervalArithmetic

function cpre_sw(Ω::Vector{IntervalBox{M,T}},
              Ω_c::Vector{IntervalBox{M,T}},
              X::IntervalBox{M,T},
              f::Vector{Function},
              ϵ::T) where {M,T<:Real}

    S = Vector{IntervalBox{M,T}}()
    N = Vector{IntervalBox{M,T}}()
    E = Vector{IntervalBox{M,T}}()
    L = Vector{IntervalBox{M,T}}(Ω)

    i = 0;
    i_s = i_n = i_e = 0
    while length(L) != 0
        x = pop!(L)

        images = [f_p(x) for f_p in f]
        if any([(!intersects_sw(image, Ω_c) && (image ⊂ X )) for image in images])
            push!(S, x)
        elseif all([!intersects_sw(image, Ω) for image in images])
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

    return (S, N, E)
end

function I∞_sw(Ω::Vector{IntervalBox{M,T}},
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
    while !isempty(N) || !isempty(E)
        Y = Ω
        Y_c = merge_best_effort([Y_c; N; E])
        (S, N, E) = cpre_sw(Y, Y_c, X, f, ϵ)
        Ω = S
        i += 1
    end

    return (S, Y_c)

end

function intersects_sw(A::IntervalBox{M,T},
                    B::Vector{IntervalBox{M,T}}) where {M,T<:Real}

    for B_i in B
        if !isempty(A ∩ B_i)
            return true
        end
    end

    return false
end


function Base.:∩(A::IntervalBox{M,T},
                 B::Vector{IntervalBox{M,T}}) where {M,T<:Real}

    res = [A ∩ B_i for B_i in B]
    res = [a for a in res if !isempty(a)]

end

function contained(A::IntervalBox{M,T},
                   B::Vector{IntervalBox{M,T}}) where {M,T<:Real}


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

    I∞(Ω, Ω_c, f, ϵ)
end

function example_pendulum_1()
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

    I∞(Ω, Ω_c, f, ϵ)
end

function example_pendulum_switched()
    n = 2
    dt = 0.01

    m = 0.2
    g = 9.8
    l = 0.3
    J = 0.006
    b = 0.1

    X = IntervalBox(-0.05..0.05, -0.01..0.01)

    U = range(-0.1,stop=0.1,step=0.02)

    f = Vector{Function}()
    for u in U
        tmp(x) = IntervalBox([x[1] + dt*x[2]; (1-dt*b/J)*x[2] + dt*(m*g*l/J*sin(x[1]) + l/J*cos(x[1])*u)])
        push!(f, tmp)
    end

    Ω = [X]
    Ω_c = [emptyinterval() × emptyinterval()]
    ϵ = 0.001

    I∞_sw(Ω, Ω_c, f, ϵ)
end

function mysqr(r::Interval{T}) where {T<:Real}
    a = r.lo;
    b = r.hi;

    af = (a >= 0);
    bf = (b >= 0);

    c = d = 0

    if (af && bf)
        c = a^2;
        d = b^2;
    elseif (af && !bf)
        c = max(a,-b)^2;
        d = 0;
    elseif (!af && bf)
        c = 0;
        d = max(-a,b)^2;
    else
        c = b^2;
        d = a^2;
    end

    return Interval(c, d)

end

function Base.:+(A::IntervalBox{M,T}, b::Vector{T}) where {M,T<:Real}

    v = Vector{Interval{T}}(undef, M)
    for i in 1:M
        v[i] = A[i] + b[i]
    end

    return IntervalBox(v)

end
