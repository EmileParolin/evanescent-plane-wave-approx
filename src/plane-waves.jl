"""
    gpw(ζ, φ; k=1)
    
Closure that defines a generalized plane wave in polar coordinates.

The closure captures the values of the evanescence parameter ``ζ``, the
angle ``φ`` and the wavenumber ``k``.
It does not include any normalization.
The function can then be evaluated at any ``(r,θ)`` point.

The generalized plane wave in polar coordinates is defined by
```math
ϕ := (r,θ) ↦ e^{ı k \\mathbf{d} ⋅ \\mathbf{x}},
```
where
```math
\\mathbf{d} = (\\cos[φ+ıζ], \\sin[φ+ıζ]),
\\qquad\\text{and}\\qquad
\\mathbf{x} = (r \\cos θ, r\\sin θ).
```
"""
function gpw(ζ,φ;k=1)
    d = [cos(φ+im*ζ), sin(φ+im*ζ)] # Direction
    return (r,θ) -> exp(im * k * [r*cos(θ), r*sin(θ)] ⋅ d)
end

"""
    approximation_set(Y, W; k=1)
    
Generate a set of generalized plane waves.

The parameters `Y` and `W` are respectively the sets of complex angles
``(ζ_j,φ_j)`` and associated weights ``ω_j`` corresponding to the generalized
plane waves, see [`gpw`](@ref)
```math
ϕ_j := (r,θ) ↦ e^{ı k \\mathbf{d}_j ⋅ \\mathbf{x}},
```
where
```math
\\mathbf{d}_j = (\\cos[φ_j+ıζ_j], \\sin[φ_j+ıζ_j]),
\\qquad\\text{and}\\qquad
\\mathbf{x} = (r \\cos θ, r\\sin θ).
```
Approximations are constructed in the form of
```math
(r, θ) ↦ \\sum_j ξ_j ω_j ϕ_j(r, θ),
```
where ``(ξ_j)_j`` is the set of unknown coefficients.

The value of the wavenumber `k` defaults to ``1`` and can be provided as an
optional parameter.
"""
function approximation_set(Y, W; k=1)
    @assert length(Y) == length(W)
    return [(r, θ) -> w * gpw(y...;k=k)(r, θ) for (y, w) in zip(Y, W)]
end

"""
    approximation_set(N; k=1)
    
Generate a set of ``N`` propagative plane waves with equispaced angles.

The equispaced angles are defined by
```math
θ_n := 2π s / R, \\qquad n = 0,…,N-1.
```

The value of the wavenumber `k` defaults to ``1`` and can be provided as an
optional parameter.
"""
function approximation_set(N; k=1)
    Y = [(0, 2π * n / N) for n in 0:N-1]
    W = ones(Int64, N) ./ sqrt(N)
    return approximation_set(Y, W; k=k)
end

"""
    approximation_set(N, qs, smpl_type; k=1)
    approximation_set(N, Q::Integer, smpl_type; k=1)
    
Generate a set of ``N`` evanescent plane waves according to some sampling type.

The second parameter can be an integer ``Q`` which controls the truncation
parameter of the [`TruncKernel`](@ref), or a set of integers `qs` to specify
only some specific terms in the expansion.
The parameter `smpl_type` should be a function name: see
[`uniform_sampling`](@ref), [`sobol_sampling`](@ref) and
[`random_sampling`](@ref).

The value of the wavenumber `k` defaults to ``1`` and can be provided as an
optional parameter.
"""
function approximation_set(N, qs, smpl_type; k=1)
    K = TruncKernel(Tuple(qs), logwz(;k=k));
    smpl = Sampler(K)
    Y, W = smpl_type(smpl, N)
    W ./= sqrt(N)
    return approximation_set(Y, W; k=k)
end
function approximation_set(N, Q::Int64, smpl_type; k=1)
    return approximation_set(N, -Q:Q, smpl_type; k=k)
end