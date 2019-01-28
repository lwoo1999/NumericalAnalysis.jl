module Eigen

using LinearAlgebra
using SparseArrays

export jacobi!, jacobi, eigenValues, eigenVectors

nonDiag(A) = [A[i, j] for i in axes(A, 1) for j in axes(A, 2) if i != j]

function findMaxNonDiag(A::Matrix{<:Real})::Tuple{Int, Int}
    elm = -Inf
    ind = (0, 0)
    for q in axes(A, 1)
        for p in 1:(q-1)
            if abs(A[p, q]) > elm
                elm = abs(A[p, q])
                ind = (p, q)
            end
        end
    end
    ind
end

function jacobi!(A::Matrix{<:Real}, e::Real)
    Q = I
    while sum(nonDiag(A).^2) > e
        p, q = findMaxNonDiag(A)
        s = (A[q,q] - A[p,p]) / 2 / A[p,q]
        t = if s == 0
            1
        else
            t1 = -s - sqrt(s^2 + 1)
            t2 = -s + sqrt(s^2 + 1)
            if abs(t1) > abs(t2)
                t2
            else
                t1
            end
        end
        c = 1 / sqrt(1 + t^2)
        d = t * c
        S = spzeros(size(A)...) + I
        S[[p, q], [p, q]] = [c d; -d c]
        Q = Q * S
        A = S' * A * S
    end
    res = []
    for i in axes(A, 1)
        push!(res, (A[i,i], collect(Q[:,i])))
    end
    res
end

function jacobi(A::Matrix{<:Real}, e::Real)::Vector{Tuple{<:Real, Vector{<:Real}}}
    A_ = copy(A)
    jacobi!(A_, e)
end

function eigenValues(A::Matrix{<:Real}, e::Real)::Vector{<:Real}
    map(first, jacobi(A, e))
end

function eigenVectors(A::Matrix{<:Real}, e::Real)::Vector{Vector{<:Real}}
    map(last, jacobi(A, e))    
end

end