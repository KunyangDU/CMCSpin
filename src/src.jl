using LinearAlgebra,LatticeUtilities,Statistics,Random
using JLD2,TimerOutputs
using CairoMakie,LaTeXStrings,ColorSchemes
# include("plot.jl")
# include("bit operations.jl")
# include("sampling.jl")

include("analysis/binning.jl")
include("dynamics/rk4.jl")
include("main/thermalize.jl")
include("main/update.jl")
include("Observables/calObs.jl")

include("lattice/abstract type.jl")
include("lattice/square.jl")
include("lattice/triangular.jl")
include("lattice/honeycomb.jl")

include("lattice/tools.jl")

include("lattice/operations.jl")
include("lattice/plot.jl")
include("lattice/geometry.jl")

include("tools.jl")

