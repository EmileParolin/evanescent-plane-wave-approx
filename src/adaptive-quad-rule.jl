"""
    function ϵsupport(f, ϵ; δ=1e-3)

Computes the ϵ-support of a function.

This is useful to compute quadratures of compactly ϵ-supported functions on
unbounded domains.
"""
function ϵsupport(f, ϵ; δ=1e-3)
    a = 0; fa = Inf; while fa > ϵ fa = abs(f(a)); a -= δ end
    b = 0; fb = Inf; while fb > ϵ fb = abs(f(b)); b += δ end
    return (a, b)
end


"""
    adaptive_G_quad(f; quadrule=gausslegendre, a=-1, b=1)
    
Computes the integral of ``f`` on ``[a,b]`` using adaptive quadrature rule.

Here adaptivity is only understood with respect to the number of nodes.
There is no-subdivision of the interval.

!!! note "Note!"
    The default implementation relies on the
    [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
    package through the `gausslegendre` quadrature rule function.
"""
function adaptive_G_quad(f; quadrule=gausslegendre, a=-1, b=1,
                         tol=1e-12, Qmin=10, Qstep=1, Qmax=10^3)
    val = 0
    err = 1
    q = Qmin
    # Estimate
    x, w = quadrule(q)
    val = sum((b-a)/2 .* w .* f.((b-a)/2 .* x .+ (a+b)/2))
    while err > tol && q < Qmax
        q += Qstep
        # Reference
        x, w = quadrule(q)
        ref = sum((b-a)/2 .* w .* f.((b-a)/2 .* x .+ (a+b)/2))
        # Error
        err = abs(val - ref) / ref
        # Update
        val = ref
    end
    if err > tol @warn "Tolerance not reached: error = $(err)." end
    if q == Qmax @warn "Maximum number $(Qmax) of nodes reached." end
    return val, err, q
end