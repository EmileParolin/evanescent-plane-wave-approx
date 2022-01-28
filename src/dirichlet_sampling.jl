"""
    samples_from_nodes(f, X)
    samples_from_nodes(f::Vector, X)
    
Evaluate `f` at the sampling nodes `X`.

If `f` is a vector, a matrix is computed.
"""
samples_from_nodes(f, X) = [f(x...) for x in X]
function samples_from_nodes(F::Vector{T}, X) where T
    return [f(x...) for x in X, f in F]
end

"""
    boundary_sampling_nodes(S)
    
The sampling nodes on the boundary of the unit disk.

The samplings nodes are then ``(1, θ_s)`` with
```math
θ_s := 2π s / S, \\qquad s = 0,…,S-1.
```
"""
boundary_sampling_nodes(S) = [(1, 2π * s / S) for s in 0:S-1]

"""
    number_of_boundary_sampling_nodes(M; η=2, P=0)
    
Number of sampling nodes necessary for an approximation set of dimension ``M``.
    
The number of sampling points on the boundary of the unit disk is given by
``S := \\max(ηM, 2P+1)``.
The (overdetermined) linear system that will be solved is then of size ``S`` by
``M``.

The oversampling parameter ``η`` defaults to ``2`` but can be provided by the
user as an optional parameter `η`.
The mode number ``P`` defaults to ``0`` but can be provided by the user as an
optional parameter `P`.
"""
number_of_boundary_sampling_nodes(M; η=2, P=0) = max(η * M, 2P+1)

"""
    Dirichlet_sampling(k, N, U; smpl_type=nothing, η=2, ϵ=1e-8, Q=(length(U)-1)//2)

Example of reconstruction of a solution surrogate via Dirichlet sampling.

The solution surrogate is defined by its vector of coefficients `U`.
The approximation set of dimension `N` can be composed of propagative plane
waves (default) or evanescent plane waves (depending on the value of
`smpl_type` and `Q` the maximum mode number assumed to be in the target).
The number of sampling points on the boundary can be controlled via the
oversampling ration `η`.
The amount of regularization in the SVD can be controlled by the regularization
parameter `ϵ`.
"""
function Dirichlet_sampling(k, U, N; smpl_type=nothing, η=2, ϵ=1e-8,
                            Q=Integer((length(U)-1)//2))
    # Maximum mode number in approximation target
    P = Integer((length(U) - 1) // 2)
    # Approximation set
    if smpl_type == nothing # Propagative plane waves
        Φ = approximation_set(N; k=k)
    else
        Φ = approximation_set(N, Q, smpl_type; k=k)
    end
    # Sampling nodes
    S = number_of_boundary_sampling_nodes(N; η=η, P=P)
    X = boundary_sampling_nodes(S)
    # Matrix and its factorization
    A = samples_from_nodes(Φ, X)
    iA = RegularizedSVDPseudoInverse(A; ϵ=ϵ)
    # Right-hand-side
    u = solution_surrogate(U; k=k)
    b = samples_from_nodes(u, X)
    # Solution
    ξ = solve_via_regularizedSVD(iA, b)
    ũ = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξ, Φ)])
    # Errors and stability estimate
    res = norm(A * ξ - b) / norm(b) # Relative residual
    nrm = norm(ξ) / norm(U)         # Stability measure
    e = (r,θ) -> ũ(r,θ) - u(r,θ)    # Absolute error function
    return u, ξ, ũ, res, nrm, e
end