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
k = 10;
nothing #hide
````

We need to define the maximum Fourier mode number in the approximation target

````@example example
P = 15;
nothing #hide
````

Next we construct the vector of coefficients in the basis `b_p`
for `p` in `[-P,P]`.

````@example example
U = zeros(ComplexF64, 2P+1);
U[P+1]     = 0.5;  # This is the constant mode (`p=0`)
U[P+1 + P] = 1im;  # This is mode `p = P`
nothing #hide
````

The target of the approximation problem can then be constructed as

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

````@example example
N = 40;
nothing #hide
````

Sampling nodes

````@example example
S = number_of_boundary_sampling_nodes(N; η=2, P=P);
X = boundary_sampling_nodes(S);
nothing #hide
````

Right-hand-side

````@example example
b = samples_from_nodes(u, X);
nothing #hide
````

## Approximation with **propagative** plane waves

### Approximation set

````@example example
Φppw = approximation_set(N; k=k);
nothing #hide
````

Matrix and its factorization

````@example example
Appw = samples_from_nodes(Φppw, X);
nothing #hide
````

### Solution

````@example example
iAppw = RegularizedSVDPseudoInverse(Appw; ϵ=1e-12);
nothing #hide
````

The coefficients are computed using the

````@example example
ξppw = solve_via_regularizedSVD(iAppw, b);
ũppw = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξppw, Φppw)]);
nothing #hide
````

### Errors and stability estimate

Relative residual

````@example example
resppw = norm(Appw * ξppw - b) / norm(b)
````

Stability measure

````@example example
nrmppw = norm(ξppw) / norm(U)
````

Absolute error function

````@example example
eppw = (r,θ) -> ũppw(r,θ) - u(r,θ);
nothing #hide
````

## Approximation with **evanescent** plane waves

### Approximation set

````@example example
Φepw = approximation_set(N, P, sobol_sampling; k=k);
nothing #hide
````

Matrix and its factorization

````@example example
Aepw = samples_from_nodes(Φepw, X);
nothing #hide
````

### Solution

````@example example
iAepw = RegularizedSVDPseudoInverse(Aepw; ϵ=1e-12);
nothing #hide
````

The coefficients are computed using the

````@example example
ξepw = solve_via_regularizedSVD(iAepw, b);
ũepw = (r,θ) -> sum([ξi * ϕi(r,θ) for (ξi, ϕi) in zip(ξepw, Φepw)]);
nothing #hide
````

### Errors and stability estimate

Relative residual

````@example example
resepw = norm(Aepw * ξepw - b) / norm(b)
````

Stability measure

````@example example
nrmepw = norm(ξepw) / norm(U)
````

Absolute error function

````@example example
eepw = (r,θ) -> ũepw(r,θ) - u(r,θ);
nothing #hide
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

