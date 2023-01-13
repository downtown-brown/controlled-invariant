import Base.merge, Base.intersect
using IntervalArithmetic

struct IntervalParameter
    intervals::Vector{Interval}
    name::String
end

struct SiviaResults{T<:IntervalBox}
    tag::String
    boxes::Vector{T}
    parameters::Vector{IntervalParameter}

    function SiviaResults(tag::AbstractString, boxes::Vector{T}, names::Vector{S}) where {T<:IntervalBox, S<:AbstractString}
        lines = [merge(project(boxes, d)) for d in 1:length(names)]
        lines = [map(x -> x.v[1], line) for line in lines]

        parameters = map(IntervalParameter, lines, names)
        new{T}(tag, boxes, parameters)
    end
end

function Base.show(io::IO, m::SiviaResults{T}) where T<:IntervalBox
    n = length(m.boxes)
    d = length(m.parameters)
    println(io, "Sivia results for $(m.tag), $n boxes, of $d parameters.")
end

function lt(A::IntervalBox{M, T}, B::IntervalBox{M, T}) where {M, T<:Real}
    for d in 1:M
        if A[d] != B[d] && A[d].lo > B[d].lo
            return false
        end
    end
    return true
end

function lt_d(A::IntervalBox{M, T}, B::IntervalBox{M, T}, d::Integer) where {M, T<:Real}
    D = 1:M
    D[1], D[d] = D[d], D[1]
    for i in D
        if A[i] != B[i] && A[i].lo > B[i].lo
            return false
        end
    end
    return true
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

function merge_pass!(A::Vector{T}) where {T<:IntervalBox}
    i = 1
    while i < length(A) - 1
        j = i + 1
        combo = merge(A[i], A[j])
        if isa(combo, T)
            A[i] = combo
            deleteat!(A, j)
            sort!(A, lt=SIVIA.lt)
        else
            i = i + 1
        end
    end
end

function merge!(A::Vector{T}) where T<:IntervalBox
    while true
        N = length(A)
        sort!(A, lt=SIVIA.lt)
        merge_pass!(A)
        length(A) == N && return
    end
end

function intersect(A::Vector{T}, B::Vector{T}, d::Integer...) where T<:IntervalBox
    Aⁱ = merge(project(A, d...))
    Bⁱ = merge(project(B, d...))

    regions = Vector{eltype(Aⁱ)}()
    for a in Aⁱ
        for b in Bⁱ
            cap = a ∩ b
            if !isempty(cap)
                push!(regions, cap)
            end
        end
    end

    return merge(regions)
end

function intersect(A::Vector{IntervalBox{M, T}}, B::Vector{IntervalBox{M, T}}) where {M, T<:Real}
    return intersect(A, B, 1:M...)
end

function project(A::Vector{T}, d...) where T<:IntervalBox
    projector(x) = IntervalBox([x.v[i] for i in d])
    projection = Vector{IntervalBox{length(d), Float64}}()
    if length(A) > 0
        append!(projection, map(projector, A))
    end

    projection
end
