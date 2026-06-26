# Continuum.jl — 连续介质力学主模块

module Continuum

include("types.jl")
include("elasticity.jl")
include("waves.jl")
include("beams.jl")
include("plates.jl")
include("failure.jl")
include("plasticity.jl")
include("fracture.jl")
include("composites.jl")
include("viscoelasticity.jl")
include("contact.jl")
include("thermoelastic.jl")
include("poroelasticity.jl")
include("conservation.jl")
include("creep.jl")
include("stability.jl")
include("micromechanics.jl")

end # module Continuum