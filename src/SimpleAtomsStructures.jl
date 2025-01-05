module SimpleAtomsStructures

# Write your package code here.
using ArgCheck
using AtomsBase
using StaticArrays
using Unitful

export AtomicPropertySystem
export CellSystem
export GeneralSystem
export GenericSystem
export SimpleAtom
export SimpleSystem
export SimpleVelocitySystem

include("abstract_systems.jl")
include("simple_atom.jl")
include("structures.jl")
include("atom_property_system.jl")
include("cell_systems.jl")
include("general_system.jl")

end
