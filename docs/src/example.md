```@meta
EditURL = "<unknown>/example/example.jl"
```

````@example example
push!(LOAD_PATH,joinpath(@__DIR__, "../src/"))
using LinearAlgebra
using StableApproxEPW
````

## Target of the approximation problem

We consider the Helmholtz solution in the unit disk, with wavenumber

````@example example
k = 5;
nothing #hide
````

We need to define the maximum Fourier mode number `P` in the approximation
target:

````@example example
P = 25;
nothing #hide
````

Next we construct the vector of coefficients in the basis `b_p`
for `p` in `[-P,P]`:

````@example example
U = zeros(ComplexF64, 2P+1);
U[P+1]     = 0.5;  # This is the constant mode (`p=0`)
U[P+1 + P] = 1im;  # This is mode `p = P`
nothing #hide
````

The target of the approximation problem can then be constructed as:

````@example example
u = solution_surrogate(U; k=k);
nothing #hide
````

``u`` can be evaluated at any `(r,θ)` point.
Alternatively, the target ``u`` could have been a single mode.
For instance, to get the circular wave with mode number `p=15`, simply set
``u = bp(15; k=k)``.
Of course, any other function (defined in polar coordinates) can be defined
by the user as target of the approximation problem.

## Reconstruction method

The approximation will be reconstructed by sampling the target on the
boundary of the unit disk.
To do so, we need to know now the dimension `N` of the approximation sets
that we are going to use:

````@example example
N = 100;
nothing #hide
````

We can determine the number of sampling nodes necessary for a successful
reconstruction, based on `N`, `P` (to avoid aliasing) and the oversampling
ratio `η`:

````@example example
S = number_of_boundary_sampling_nodes(N; η=2, P=P);
nothing #hide
````

The (equispaced) boundary nodes are then constructed as:

````@example example
X = boundary_sampling_nodes(S);
nothing #hide
````

The right-hand-side of the linear system can then be readily constructed:

````@example example
b = samples_from_nodes(u, X);
nothing #hide
````

## Approximation with **propagative** plane waves

We construct the approximation set of PPW:

````@example example
Φppw = approximation_set(N; k=k);
nothing #hide
````

Matrix and its factorization

````@example example
Appw = samples_from_nodes(Φppw, X);
iAppw = RegularizedSVDPseudoInverse(Appw; ϵ=1e-14);
nothing #hide
````

The coefficients of the approximation are computed:

````@example example
ξppw = solve_via_regularizedSVD(iAppw, b);
ũppw = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξppw, Φppw)]);
nothing #hide
````

Absolute error function

````@example example
eppw = (r,θ) -> ũppw(r,θ) - u(r,θ); eppw(1,π/2)
````

Let's compute the relative residual:

````@example example
resppw = norm(Appw * ξppw - b) / norm(b)
````

We did not obtained any accuracy whatsoever!
The reason is that the coefficients are too large:

````@example example
nrmppw = norm(ξppw) / norm(U)
````

## Approximation with **evanescent** plane waves

We construct the approximation set of EPW:

````@example example
Φepw = approximation_set(N, P, sobol_sampling; k=k);
nothing #hide
````

Matrix and its factorization

````@example example
Aepw = samples_from_nodes(Φepw, X);
iAepw = RegularizedSVDPseudoInverse(Aepw; ϵ=1e-14);
nothing #hide
````

The coefficients of the approximation are computed:

````@example example
ξepw = solve_via_regularizedSVD(iAepw, b);
ũepw = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξepw, Φepw)]);
nothing #hide
````

Absolute error function

````@example example
eepw = (r,θ) -> ũepw(r,θ) - u(r,θ); eepw(1,π/2)
````

Let's compute the relative residual:

````@example example
resepw = norm(Aepw * ξepw - b) / norm(b)
````

We get almost 8 digits of accuracy!
The size of the coefficients remains quite high:

````@example example
nrmepw = norm(ξepw) / norm(U)
````

### Down to machine precision

Let's double the number of EPW

````@example example
N = 200
````

The above approximation process can be obtained much more rapidly with the
following convenience function

````@example example
_, ξepw, ũepw, resepw, nrmepw, eepw = Dirichlet_sampling(k, U, N; smpl_type=sobol_sampling);
nothing #hide
````

One can check that we obtain now a relative residual very close to machine
precision (13 digits of accuracy):

````@example example
resepw
````

The size of the coefficients is also greatly reduced:

````@example example
nrmepw
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

