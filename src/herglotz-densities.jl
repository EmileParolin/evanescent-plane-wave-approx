"""
    ãp(p)
    
Closure for the Herglotz polynomial.

The closure captures the values of the mode number ``p`` and does not include any
normalization.
The function can then be evaluated at any ``(ζ,φ)`` point in the complex strip.
Here ``ζ`` is the evanescence parameter and ``φ`` is the angle indicating the
direction of propagation.

The Herglotz polynomials are defined by
```math
ã_p := (ζ,φ) ↦ e^{p (ζ + ıφ)}, \\qquad ∀ p ∈ \\mathbb{Z}.
```
"""
ãp(p) = (ζ, φ) -> exp(p * (ζ + im * φ))
logãp(p) = (ζ, φ) -> p * (ζ + im * φ)

"""
    wz(; k=1, z=1/4)
    
Closure to define the weight function.

The closure captures the values of the wavenumber ``k`` and a parameter ``z``.
The function can then be evaluated at any ``ζ`` point.

The weight function is defined as
```math
w := ζ ↦ e^{z |ζ| - k \\sinh|ζ|}, \\qquad ∀ ζ ∈ \\mathbb{R}.
```
"""
wz(;k=1, z=1/4) = ζ -> exp(z * abs(ζ) - k * sinh(abs(ζ)))
logwz(;k=1, z=1/4) = ζ -> z * abs(ζ) - k * sinh(abs(ζ))

"""
    αp(p, logweight)
    
Computation of the normalization constant of the Herglotz polynomial.

The Herglotz polynomial is defined by ``ã_p``, see the function [`ãp`](@ref).
The normalization constant is defined by
```math
α_p^{-2} := 2π \\int_{\\mathbb{R}} |ã_p|^{2} w^2 \\mathrm{d}\\zeta,
```
where ``w`` is the weight function.
The weight function ``w`` can be taken to be the function [`wz`](@ref).

The implementation relies on the convenience functions
[`ϵsupport`](@ref) and [`adaptive_G_quad`](@ref).

!!! danger "Beware!"
    The second argument expects the ``\\log`` of the weight function ``w``.
    This is for numerical stability reasons.

!!! warning "Warning!"
    The implementation uses some ad-hoc adaptive Gauss quadrature routine.
    There is probably some room for improvement in the implementation.
"""
function αp(p, logweight; ϵ=1e-14)
    # Definition of the integrand
    intgd = ζ -> abs(exp(2 * logãp(p)(ζ, 0) + 2 * logweight(ζ)))
    # ϵ-support of the integrand
    ζm∞, ζp∞ = ϵsupport(intgd, ϵ)
    # Adaptive Gauss quadrature on the ϵ-support
    αmquad = adaptive_G_quad(intgd; a=ζm∞, b=0)
    αpquad = adaptive_G_quad(intgd; a=0, b=ζp∞)
    # Sanity checks
    if !(αmquad[1] < Inf || αpquad[1] < Inf)
        msg = "p = $(p) too large: a_p not integrable:"
        msg *= "\n$(ζm∞)<x<0 (val, err, q) = $(αmquad)"
        msg *= "\n0<x<$(ζp∞) (val, err, q) = $(αpquad)"
        error(msg)
    end
    if (αmquad[2] + αpquad[2]) > 1e-8
        msg = "p = $(p): normalization constant may be inaccurate:"
        msg *= "\n$(ζm∞)<x<0 (val, err, q) = $(αmquad)"
        msg *= "\n0<x<$(ζp∞) (val, err, q) = $(αpquad)"
        @warn msg
    end
    # Normalization using that the integrand is constant w.r.t φ ∈ [0,2π]
    α = 1 / sqrt(2π * (αmquad[1] + αpquad[1]))
    return α
end

"""
    ap(p, logweight)
    
Closure for the normalized Herglotz polynomial.

The closure captures the values of the mode number ``p``.
The Herglotz polynomial is defined by ``ã_p``, see the function [`ãp`](@ref).
The normalization is given by the second argument, as defined in the
documentation of the function [`αp`](@ref) and is precomputed once for all.
The function can then be evaluated at any ``(ζ,φ)`` point in the complex strip.
Here ``ζ`` is the evanescence parameter and ``φ`` is the angle indicating the
direction of propagation.

The normalized Herglotz polynomials are defined by
```math
a_p := α_p ã_p, \\qquad ∀ p ∈ \\mathbb{Z}.
```
"""
function ap(p, logweight)
    α = αp(p, logweight)
    return (ζ, φ) -> α * ãp(p)(ζ, φ)
end