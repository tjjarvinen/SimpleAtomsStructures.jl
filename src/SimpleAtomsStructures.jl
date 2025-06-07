module SimpleAtomsStructures

using ArgCheck
using AtomsBase
using Compat
using LinearAlgebra: cross, dot, norm, I
using Rotations
using StaticArrays
using Unitful
using Reexport

# Structures
export AtomicPropertySystem
export AtomicPropertySystemView
export CellSystem
export CellSystemView
export GenericSystem
export SimpleAtom
export SimpleSystem
export SimpleSystemView
export SimpleVelocitySystem
export SimpleVelocitySystemView
export SystemView
export VelocityTrajectory

# utilities
export add_systems
export bond_angle
export cell_matrix
export center_of_mass
export dihedral_angle
export distance
export distance_vector
export fractional_coordinates
export fractional_coordinates_as_matrix
export inv_cell
export position_as_matrix
export rotate_system
export rotate_system!
export system_view
export translate_system
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
include("subsystem_view.jl")
include("utils.jl")

#submodule
include("trajectory/AtomsTrajectories.jl")

@reexport using .AtomsTrajectories

end
