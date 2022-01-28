"""
    TruncKernel(P::Integer, logweight)
    TruncKernel(ps::NTuple{N,Int64}, logweight) where N

Structure to hold precomputed normalization constants of Herglotz polynomials.

This is a convenience structure to hold precomputed normalization constants
[`αp`](@ref) of Herglotz polynomials which are expensive to compute.

All the normalized Herglotz polynomials ``a_p`` (see [`ap`](@ref)) for ``p`` in
the argument `ps` (noted ``\\mathcal{P}``) are precomputed and stored in the
attribute `as` which is a `Tuple`.
The `logweight` function is used to compute the normalization constants.
If only an `Integer` ``P`` is given as first argument, the terms in the kernel
are the ``a_p`` functions for ``p`` ranging from ``-P`` to ``P``.

A kernel can be evaluated at `(x,y)` where `x=(ζx,φx)` and `y=(ζy,φy)` are
points in the complex strip.

```math
(\\mathbf{x},\\mathbf{y}) ↦ K(\\mathbf{x}, \\mathbf{y})
= \\sum_{p \\in \\mathcal{P}} \\overline{a_{p}(\\mathbf{x})} a_{p}(\\mathbf{y}).
```
"""
struct TruncKernel{N}
    ps::NTuple{N,Int64}
    logweight
    as::NTuple{N} # Tuple containing the a_p functions for p in ps
end
function TruncKernel(ps::NTuple{N,Int64}, logweight) where N
    # Computing the weights
    αs = [αp(p, logweight) for p in ps]
    # ãp's Herglotz polynomials
    ãs = [ãp(p) for p in ps]
    # Normalized ap's Herglotz polynomials
    as = [(ζ,φ) -> α * ã(ζ,φ) for (α,ã) in zip(αs, ãs)]
    return TruncKernel(ps, logweight, Tuple(as))
end
TruncKernel(P::Int64, logweight) = TruncKernel(Tuple(-P:P), logweight)
(K::TruncKernel)(x, y) = sum([ajoint(a(x...))*a(y...) for a in K.as])

"""
    kernel_diagonal(K::TruncKernel)

Closure that defines the evaluation function of the diagonal of a
[`TruncKernel`](@ref).

```math
\\mathbf{y} ↦ K(\\mathbf{y}, \\mathbf{y})
= \\sum_{p \\in \\mathcal{P}} |a_{p}(\\mathbf{y})|^2,
\\qquad\\mathbf{y}=(ζ,φ).
```
"""
function kernel_diagonal(K::TruncKernel)
    return (ζ, φ) -> sum([abs(a(ζ, φ))^2 for a in K.as])
end

"""
    probabilitydensityfunction(K::TruncKernel)

Probability density function defined by a [`TruncKernel`](@ref).

```math
(ζ, φ) ↦ \\frac{1}{|\\mathcal{P}|} K(\\mathbf{y}, \\mathbf{y})
= \\frac{1}{|\\mathcal{P}|} \\sum_{p \\in \\mathcal{P}} |a_{p}(\\mathbf{y})|^2,
\\qquad\\mathbf{y}=(ζ,φ).
```

Notice that it is indeed a PDF since all the ``a_p`` functions are orthonormal.
"""
function probabilitydensityfunction(K::TruncKernel)
    return (ζ, φ) -> kernel_diagonal(K)(ζ, φ) / length(K.ps)
end
function probabilitydensityfunction_weighted(K::TruncKernel)
    return (ζ, φ) -> probabilitydensityfunction(K)(ζ, φ) * exp(2 * K.logweight(ζ))
end

"""
    cumulativedensityfunction(K::TruncKernel)

Cumulative density function defined by a [`TruncKernel`](@ref).

The probability density function is provided by the function
[`probabilitydensityfunction`](@ref).

```math
ζ ↦ \\frac{2π}{|\\mathcal{P}|} \\int_{-\\infty}^{ζ} K(ζ,ζ) \\mathrm{d}ζ
= \\frac{2π}{|\\mathcal{P}|} \\int_{-\\infty}^{ζ} \\sum_{p \\in \\mathcal{P}} |a_{p}|^2(ζ) \\mathrm{d}ζ,
```
where, abusing notations, we used ``K(ζ,ζ) = K(\\mathbf{y},\\mathbf{y})``
and ``|a_{p}|^2(ζ) = |a_{p}|^2(\\mathbf{y})`` for any ``\\mathbf{y} = (ζ,φ)``
since they both are quantities independent of ``φ``.

The implementation relies on the convenience functions
[`ϵsupport`](@ref) and [`adaptive_G_quad`](@ref).

!!! danger "Beware!"
    The implementation *assumes* that the probability density function is
    constant with respect to the second variable ``φ``.
    This is the reason why the resulting CDF is a univariate function.
"""
function cumulativedensityfunction(K::TruncKernel; tol=1e-12)
    pdf = ζ -> probabilitydensityfunction_weighted(K)(ζ, 0)
    # Chunking the support of the integrand for some h-adaptivity
    ϵs = pdf(0) .* 10. .^ (-12:3:0)
    supps = sort(vcat(0, [collect(ϵsupport(pdf, ϵ)) for ϵ in ϵs]...))
    as = supps[1:end-1]
    bs = supps[2:end]
    # CDF evaluation on whole real line and creation of quadrature rule
    cdf∞ = Float64(0)
    xs = Vector{Float64}[]; ws = Vector{Float64}[];
    for (a, b) in zip(as, bs)
        val, _, q = adaptive_G_quad(pdf; a=a, b=b, tol=tol)
        cdf∞ += 2π * val # using that CDF constant w.r.t φ ∈ [0,2π]
        x, w = gausslegendre(q)
        push!(xs, x); push!(ws, w)
    end
    # Sanity checks
    if abs(supps[1] + supps[end]) > 1e-12
        @warn "PDF ϵ-support = [$(supps[1]),$(supps[end])] not symmetric."
    end
    if abs(cdf∞ - 1) > 1e-8
        @warn "PDF not normalized ∫pdf = $(cdf∞)."
    end
    # CDF evaluation function
    function cdf(ζ)
        c = Float64(0)
        for (a, b, x, w) in zip(as, bs, xs, ws)
            if a < ζ
                b = min(b, ζ)
                c += 2π * sum((b-a)/2 .* w .* pdf.((b-a)/2 .* x .+ (a+b)/2))
            end
        end
        return c
    end
    return cdf
end

"""
    Sampler(K::TruncKernel)

Structure to compute sampling nodes and weights.

The nodes and weights are sampled according to the probability density function
defined by a [`TruncKernel`](@ref).

!!! note "Note!"
    The implementation uses the secant method as root finding algorithm to
    compute the pre-image of the CDF.
    There is probably some room for improvement in the implementation,
    specially since we know that the first derivative of the CDF is the PDF.
    The Newton-Raphson method fails to converge for instance.
"""
struct Sampler
    K::TruncKernel
    pdf
    cdf
    supp
end
function Sampler(K::TruncKernel)
    pdf = ζ -> probabilitydensityfunction_weighted(K)(ζ, 0)
    cdf = cumulativedensityfunction(K)
    supp = ϵsupport(pdf, pdf(0) * 10. ^ -12)
    return Sampler(K, pdf, cdf, supp)
end

"""
    inversion(smpl::Sampler, u)
    
Computes the pre-image ``ζ`` of `u` under the cumulative density function.

This function allows to perform *Inversion Transform Sampling*.

This is a univariate inversion function because the cdf is constant with
respect to the other variable ``φ``.
"""
function inversion(smpl::Sampler, u)
    f = x -> smpl.cdf(x) - u # Objective function
    ζ = begin 
        try
            find_zero((f, smpl.pdf), 0, Roots.Secant())
        catch # Fall back to slower but robust method if not converging
            find_zero((f, smpl.pdf), smpl.supp)
        end
    end
    return ζ
end

"""
    weight(smpl::Sampler, node)
    
Computes the weight associated to a node ``(ζ,φ)``.

The weight at ``\\mathbf{y}=(ζ,φ)`` is defined by
```math
ω(\\mathbf{y}) := \\sqrt{\\frac{|\\mathcal{P}|}{\\sum_{p \\in \\mathcal{P}} |a_{p}(\\mathbf{y})|^2}},
\\qquad\\mathbf{y}=(ζ,φ).
```
"""
function weight(smpl::Sampler, node)
    return 1 / sqrt(probabilitydensityfunction(smpl.K)(node...))
end

"""
    random_sampling(smpl::Sampler)
    random_sampling(smpl::Sampler, M)

Draw `M` random samples according to the distribution stored in `smpl`.
"""
function random_sampling(smpl::Sampler)
    node = (inversion(smpl, rand()), 2π * rand())
    wght = weight(smpl, node)
    return node, wght
end
function random_sampling(smpl::Sampler, M)
    xw = [random_sampling(smpl) for _ in 1:M]
    nodes = [x for (x,_) in xw]
    wghts = [w for (_,w) in xw]
    return nodes, wghts
end

"""
    uniform_sampling(smpl::Sampler, M)

Draw `M` deterministic samples according to the distribution stored in `smpl`.

This has a tensor-like structure, with the same number of samples in both
directions.
"""
function uniform_sampling(smpl::Sampler, M)
    Ms = Integer(ceil(sqrt(M)))
    φs = range(0, stop=2π, length=Ms+1)[1:end-1]
    uζs = range(0, stop=1, length=Ms+1)[1:end-1]
    uζs = collect(uζs) .+ uζs[2]/2
    nodes = [(inversion(smpl, u), φ) for u in uζs, φ in φs][:]
    wghts = weight.(Ref(smpl), nodes)
    return nodes, wghts
end

"""
    sobol_sampling(smpl::Sampler, M)

Draw `M` quasi-random samples according to the distribution stored in `smpl`.

This function uses Sobol sequences computed according to the package
[QuasiMonteCarlo.jl](https://github.com/SciML/QuasiMonteCarlo.jl).
"""
function sobol_sampling(smpl::Sampler, M)
    s = QuasiMonteCarlo.sample(M,[0,0],[1,2π],SobolSample())
    ζs = s[1,:]
    φs = s[2,:]
    nodes = [(inversion(smpl, ζ), φ) for (ζ, φ) in zip(ζs, φs)][:]
    wghts = weight.(Ref(smpl), nodes)
    return nodes, wghts
end