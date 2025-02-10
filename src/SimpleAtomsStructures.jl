module SimpleAtomsStructures

# Write your package code here.
using ArgCheck
using AtomsBase
using Compat
using LinearAlgebra: cross, dot, norm, I
using Rotations
using StaticArrays
using Unitful

export AtomicPropertySystem
export CellSystem
export GenericSystem
export SimpleAtom
export SimpleSystem
export SimpleVelocitySystem

# utilities
export angled
export cell_matrix
export center_of_mass
export dihedral_angle
export dihedral_angled
export distance
export distance_vector
export fractional_coordinates
export fractional_coordinates_as_matrix
export inv_cell
export position_as_matrix
export rotate_system!
export translate_system!
export wrap_coordinates!


@compat public AbstractCompositeSystem
@compat public AbstractIsolatedSystem
@compat public AbstractSimpleSystem
@compat public GeneralSystem

include("abstract_systems.jl")
include("simple_atom.jl")
include("simple_systems.jl")
include("atom_property_system.jl")
include("cell_systems.jl")
include("general_system.jl")
include("utils.jl")

end
