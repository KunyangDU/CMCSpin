using LinearAlgebra,LatticeUtilities,Statistics,Random
using JLD2,TimerOutputs
using CairoMakie,LaTeXStrings,ColorSchemes

include("lattice/abstract type.jl")
include("type/abstract type.jl")

include("type/control.jl")
include("type/operator.jl")
include("type/state.jl")

include("analysis/binning.jl")
include("dynamics/rk4.jl")

include("algorithm/thermalize.jl")
include("algorithm/SD.jl")
include("algorithm/SA.jl")
include("algorithm/update.jl")

include("Hamiltonian/addIntr.jl")
include("Hamiltonian/operations.jl")

include("lattice/simple lattice.jl")
include("lattice/composite lattice.jl")
include("lattice/tools.jl")
include("lattice/operations.jl")
include("lattice/plot.jl")
include("lattice/geometry.jl")

include("tools.jl")

