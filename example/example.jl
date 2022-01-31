push!(LOAD_PATH,joinpath(@__DIR__, "../src/"))
using LinearAlgebra
using StableApproxEPW

# ## Target of the approximation problem

# We consider the Helmholtz solution in the unit disk, with wavenumber
k = 10;

# We need to define the maximum Fourier mode number in the approximation target
P = 15;

# Next we construct the vector of coefficients in the basis `b_p`
# for `p` in `[-P,P]`.
U = zeros(ComplexF64, 2P+1);
U[P+1]     = 0.5;  # This is the constant mode (`p=0`)
U[P+1 + P] = 1im;  # This is mode `p = P`

# The target of the approximation problem can then be constructed as
u = solution_surrogate(U; k=k);
# ``u`` can be evaluated at any `(r,θ)` point.
# Alternatively, the target ``u`` could have been a single mode.
# For instance, to get the circular wave with mode number `p=15`, simply set
# ``u = bp(15; k=k)``.
# Of course, any other function (defined in polar coordinates) can be defined
# by the user as target of the approximation problem.

# ## Reconstruction method

N = 40;

# Sampling nodes
S = number_of_boundary_sampling_nodes(N; η=2, P=P);
X = boundary_sampling_nodes(S);

# Right-hand-side
b = samples_from_nodes(u, X);


# ## Approximation with **propagative** plane waves

# ### Approximation set

Φppw = approximation_set(N; k=k);
# Matrix and its factorization
Appw = samples_from_nodes(Φppw, X);


# ### Solution

iAppw = RegularizedSVDPseudoInverse(Appw; ϵ=1e-12);

# The coefficients are computed using the 
ξppw = solve_via_regularizedSVD(iAppw, b);
ũppw = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξppw, Φppw)]);


# ### Errors and stability estimate

# Relative residual
resppw = norm(Appw * ξppw - b) / norm(b)
# Stability measure
nrmppw = norm(ξppw) / norm(U)
# Absolute error function
eppw = (r,θ) -> ũppw(r,θ) - u(r,θ);


# ## Approximation with **evanescent** plane waves

# ### Approximation set

Φepw = approximation_set(N, P, sobol_sampling; k=k);
# Matrix and its factorization
Aepw = samples_from_nodes(Φepw, X);


# ### Solution

iAepw = RegularizedSVDPseudoInverse(Aepw; ϵ=1e-12);

# The coefficients are computed using the 
ξepw = solve_via_regularizedSVD(iAepw, b);
ũepw = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξepw, Φepw)]);


# ### Errors and stability estimate

# Relative residual
resepw = norm(Aepw * ξepw - b) / norm(b)
# Stability measure
nrmepw = norm(ξepw) / norm(U)
# Absolute error function
eepw = (r,θ) -> ũepw(r,θ) - u(r,θ);