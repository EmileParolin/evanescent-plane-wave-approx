"""
    b̃p(p; k=1)
    
Closure for the circular wave function in polar coordinates.

The closure captures the values of the mode number ``p`` as well as the
wavenumber ``k`` and does not include any normalization.
The function can then be evaluated at any ``(r,θ)`` point.

The circular waves in polar coordinates are defined by
```math
b̃_p := (r,θ) ↦ J_{p}(kr) e^{ı p θ}, \\qquad ∀ p ∈ \\mathbb{Z}.
```
"""
b̃p(p; k=1) = (r,θ) -> besselj(p, k*r) * exp(im * p * θ)

"""
    βp(p; k=1)
    
Computation of the normalization constant of the circular waves.

The circular wave is defined by [`b̃p`](@ref) and the normalization is done
using the ``k``-weighted ``H^1`` norm.

By definition
```math
    \\| b̃_p \\|_{H^{1}}^2 = \\| b̃_p \\|_{L^{2}}^2 + k^{-2} \\| \\nabla b̃_p \\|_{L^{2}}^2.
```
Using integration by parts,
```math
    \\| \\nabla b̃_p \\|_{L^{2}}^2 = k^2 \\| b̃_p \\|_{L^{2}}^2 + (\\partial_n b̃_p , b̃_p).
```
On the one hand
```math
    \\| b̃_p \\|_{L^{2}}^2 = \\pi (J_p^2(k) - J_{p-1}(k)J_{p+1}(k)).
```
On the other hand
```math
    (\\partial_n b̃_p , b̃_p) = \\pi k (J_{p-1}(k) - J_{p+1}(k))J_p(k).
```
Hence
```math
    \\| b̃_p \\|_{H^{1}}^2 = 2\\pi (J_p^2(k) - J_{p-1}(k)J_{p+1}(k))
    + \\pi k^{-1} (J_{p-1}(k) - J_{p+1}(k))J_p(k).
```
"""
function βp(p; k=1)
    L2nrm2 = π * (besselj(p,k)^2 - besselj(p-1,k) * besselj(p+1,k))
    H1nrm2 = 2 * L2nrm2 + (π/k) * (besselj(p-1,k) - besselj(p+1,k)) * besselj(p,k)
    return 1 / sqrt(H1nrm2)
end

"""
    bp(p; k=1)
    
Closure for the normalized circular waves.

The closure captures the values of the mode number ``p`` and the wavenumber
``k``.
The normalization is the default normalization defined in the function
[`βp`](@ref) and is precomputed once for all.
The function can then be evaluated at any ``(r,θ)`` point.

The normalized circular waves are defined by
```math
b_p := β_p b̃_p, \\qquad ∀ p ∈ \\mathbb{Z}.
```
"""
function bp(p; k=1)
    β = βp(p; k=k)
    return (r,θ) -> β * b̃p(p; k=k)(r,θ)
end

"""
    solution_surrogate(U; k=1)
    
Computes a solution surrogate from a set of coefficients `U`.

The parameter `U` is a vector of coefficients ``u_p`` for ``p`` ranging from
``-P`` to ``P`` (hence of size ``2P+1``).
The solution surrogate is then
```math
\\mathbf{x} ↦ \\sum_{|p| \\leq P} u_p b_{p}(\\mathbf{x}).
```
"""
function solution_surrogate(U; k=1)
    P = (length(U) - 1) // 2
    return (r, θ) -> sum([u * bp(p; k=k)(r, θ) for (p, u) in zip(-P:P, U)])
end