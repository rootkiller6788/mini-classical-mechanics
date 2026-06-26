# Chaos.jl — 混沌动力学主模块
# 零外部依赖 (仅 LinearAlgebra stdlib)
#
# 参考: Strogatz Nonlinear Dynamics, Goldstein Ch.11

module Chaos

include("types.jl")
include("flows.jl")
include("maps.jl")
include("lyapunov.jl")
include("bifurcation.jl")
include("fractals.jl")
include("synchronization.jl")
include("embedding.jl")
include("control.jl")
include("networks.jl")
include("recurrence.jl")
include("stochastic.jl")
include("timeseries.jl")

end # module Chaos
