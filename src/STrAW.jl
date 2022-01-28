module STrAW

using LinearAlgebra, SparseArrays, SpecialFunctions
using FastGaussQuadrature, Roots
using QuasiMonteCarlo

include("adaptive-quad-rule.jl")
export ϵsupport, adaptive_G_quad

include("circular-waves.jl")
export b̃p, βp, bp, solution_surrogate

include("herglotz-densities.jl")
export ãp, αp, ap, wz, logwz

include("plane-waves.jl")
export gpw, approximation_set

include("parametric_sampling.jl")
export TruncKernel, kernel_diagonal
export probabilitydensityfunction, probabilitydensityfunction_weighted
export cumulativedensityfunction
export Sampler, weight, inversion
export random_sampling, uniform_sampling, sobol_sampling

include("regularizedSVD.jl")
export RegularizedSVDPseudoInverse, condition_number, solve_via_regularizedSVD

include("dirichlet_sampling.jl")
export samples_from_nodes, boundary_sampling_nodes
export number_of_boundary_sampling_nodes
export Dirichlet_sampling

end
