# CelestialMechanics.jl — 天体力学主模块
# 零外部依赖 (仅 LinearAlgebra stdlib)
#
# 参考: Goldstein Ch.3, Murray & Dermott, Vallado

module CelestialMechanics

include("types.jl")
include("kepler.jl")
include("two_body.jl")
include("perturbations.jl")
include("three_body.jl")
include("mission.jl")

end # module CelestialMechanics
