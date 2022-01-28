"""
    RegularizedSVDPseudoInverse(A; ϵ=1e-8)

Structure able to solve linear system `A` using a regularized SVD.

The threshold `ϵ` controls the regularization. 
Let ``σ_{\\max}`` be the larger singular value.
A singular value ``σ`` such that ``σ \\leq ϵ σ_{\\max}`` is approximated by
zero when solving the least-squares problem.
"""
struct RegularizedSVDPseudoInverse
    A
    ϵ
    F
    I
    Uϵ
    iΣϵ
    Vϵ
    function RegularizedSVDPseudoInverse(A; ϵ=1e-8)
        F = svd(A)
        I = argmin(F.S .> ϵ * F.S[1]) - 1
        if I == 0 I = length(F.S) end
        iΣϵ = spdiagm(0 => 1 ./ F.S[1:I])
        return new(A, ϵ, F, I, F.U[:,1:I], iΣϵ, F.V[:,1:I]) 
    end
end

"""
    condition_number(B::RegularizedSVDPseudoInverse)
    
Computes condition number of the matrix for which `B` was constructed.

The condition number is defined as the ratio of the largest over the smallest
singular values of the matrix ``A``.
```math
κ := \\frac{σ_{\\max}}{σ_{\\min}}.
```
"""
function condition_number(B::RegularizedSVDPseudoInverse)
    return B.F.S[1] / B.F.S[end]
end

"""
    solve_via_regularizedSVD(B::RegularizedSVDPseudoInverse, b)

Solve the linear system ``Ax=b`` using the regularized SVD of ``A``.
"""
function solve_via_regularizedSVD(B::RegularizedSVDPseudoInverse, b)
    return B.Vϵ * ( B.iΣϵ * ( adjoint(B.Uϵ) * b ))
end