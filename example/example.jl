push!(LOAD_PATH,joinpath(@__DIR__, "../src/"))
using LinearAlgebra
using StableApproxEPW

# ## Target of the approximation problem

# We consider the Helmholtz solution in the unit disk, with wavenumber
k = 5;
# We need to define the maximum Fourier mode number `P` in the approximation
# target:
P = 25;
# The target approximation will be of the form
# ```math
# u = \mathbf{x} ↦ \sum_{|p| \leq P} u_p b_{p}(\mathbf{x}).
# ```
# Next we construct the vector of coefficients in the basis `b_p`
# for `p` in `[-P,P]`:
U = zeros(ComplexF64, 2P+1);
U[P+1]     = 0.5;  # This is the constant mode (`p=0`)
U[P+1 + P] = 1im;  # This is mode `p = P`
# Here the only non-zero coefficients are ``u_0 = \frac{1}{2}`` and
# ``u_P = \imath`` so that
# ```math
# u = \frac{1}{2} b_{0} + \imath b_{P}.
# ```
# The target of the approximation problem can then be constructed as:
u = solution_surrogate(U; k=k);
# ``u`` can be evaluated at any `(r,θ)` point.
# Alternatively, the target ``u`` could have been a single mode.
# For instance, to get the circular wave with mode number `p=15`, simply set
# ``u = bp(15; k=k)``.
# Of course, any other function (defined in polar coordinates) can be defined
# by the user as target of the approximation problem.

# ## Reconstruction method

# The approximation will be reconstructed by sampling the target on the
# boundary of the unit disk.
# To do so, we need to know now the dimension `N` of the approximation sets
# that we are going to use:
N = 100;
# We can determine the number of sampling nodes necessary for a successful
# reconstruction, based on `N`, `P` (to avoid aliasing) and the oversampling
# ratio `η`:
S = number_of_boundary_sampling_nodes(N; η=2, P=P);
# The (equispaced) boundary nodes are then constructed as:
X = boundary_sampling_nodes(S);
# The right-hand-side of the linear system can then be readily constructed:
b = samples_from_nodes(u, X);


# ## Approximation with **propagative** plane waves

# We construct the approximation set of PPW:
Φppw = approximation_set(N; k=k);
# The matrix and its (SVD) factorization
Appw = samples_from_nodes(Φppw, X);
iAppw = RegularizedSVDPseudoInverse(Appw; ϵ=1e-14);
# The coefficients of the approximation are computed:
ξppw = solve_via_regularizedSVD(iAppw, b);
ũppw = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξppw, Φppw)]);
# Absolute error function
eppw = (r,θ) -> ũppw(r,θ) - u(r,θ); eppw(1,π/2)
# Let's compute the relative residual:
resppw = norm(Appw * ξppw - b) / norm(b)
# We did not obtained any accuracy whatsoever!
# The reason is that the coefficients are too large:
nrmppw = norm(ξppw) / norm(U)


# ## Approximation with **evanescent** plane waves

# We construct the approximation set of EPW:
Φepw = approximation_set(N, P, sobol_sampling; k=k);
# It is possible to choose other types of sampling methods, instead of
# `sobol_sampling`, for instance `uniform_sampling` and `random_sampling`.
# The matrix and its (SVD) factorization
Aepw = samples_from_nodes(Φepw, X);
iAepw = RegularizedSVDPseudoInverse(Aepw; ϵ=1e-14);
# The coefficients of the approximation are computed:
ξepw = solve_via_regularizedSVD(iAepw, b);
ũepw = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξepw, Φepw)]);
# Absolute error function
eepw = (r,θ) -> ũepw(r,θ) - u(r,θ); eepw(1,π/2)
# Let's compute the relative residual:
resepw = norm(Aepw * ξepw - b) / norm(b)
# We get almost 8 digits of accuracy!
# The size of the coefficients remains quite high:
nrmepw = norm(ξepw) / norm(U)

# ### Down to machine precision

# Let's double the number of EPW
N = 200
# The above approximation process can be obtained much more rapidly with the
# following convenience function
_, ξepw, ũepw, resepw, nrmepw, eepw = Dirichlet_sampling(k, U, N; smpl_type=sobol_sampling);
# One can check that we obtain now a relative residual very close to machine
# precision (13 digits of accuracy):
resepw
# The size of the coefficients is also greatly reduced:
nrmepw