# SimpleAtomsStructures

[![Build Status](https://github.com/tjjarvinen/SimpleAtomsStructures.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tjjarvinen/SimpleAtomsStructures.jl/actions/workflows/CI.yml?query=branch%3Amain)


New implementations for AtomsBase system structures

## Atom Structure

New atom structure `SimpleAtom` is a bits type structure that aims to be a simple small structure.

You can create `SimpleAtom` in one of the following ways

```julia
SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å")
SimpleAtom( 1, [0.0, 0.0, 0.0]u"Å") # same as above
SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å", [0.1, 0.0, 0.0]u"Å/s"; mass = 16.0u"u", charge = -1.0u"q")
SimpleAtom(ChemicalSpecies(:H), [0.0, 0.0, 0.0]u"Å")
SimpleAtom( :O => [1.0, 0.0, 0.0]u"Å" )
```


Comparison to AtomsBase `Atom`
```julia-repl
julia> ab_atom = AtomsBase.Atom( :O, [1.0, 0.0, 0.0]u"Å" )
Atom(O, Z = 8, m = 15.999 u):
    position          : [1,0,0]u"Å"
    species           : O

julia> sa = SimpleAtom( :O, [1.0, 0.0, 0.0]u"Å" )
SimpleAtom(O, Z = 8, m = 15.999 u):
    position          : [1,0,0]u"Å"
    species           : O

julia> Base.summarysize(ab_atom)
456

julia> Base.summarysize(sa)
32

julia> isbits(ab_atom)
false

julia> isbits(sa)
true
```

## System Structures

Under the hood there are several different system structures that differ on what data they store. But you only need to call `generic_system` to build all of them.
Depending on what information you provide you get different structure

```julia
# Build based on vector of atoms
# generic_system(atoms::AbstractVector{SimpeAtom}; kwargs...)
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"), 
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]
)

# Same but added key for energy
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")];
    energy = 10.0u"eV"
)

# Add cell to the system
generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]; 
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"], 
    periodicity = (true, true, true)
)

# Create a system from an array of pairs
# generic_system(AbstractVector(<:Pair); kwargs...)
sys = generic_system([:H => [0.0, 0.0, 0.0]u"Å", :O => [1.0, 0.0, 0.0]u"Å"])

# Create a system vectors of atom symbols, positions and velocities
# geric_system(spc, pos, [vel]; kwargs...)
sys = generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"],
    [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"];
)
```

### Build System from Other Systems

You can build system from other system and modify system global features

```julia
# Form a copy of old system
# gneric_system(old_sys; kwargs...)
new_sys = generic_system(old_sys)

# Copy system and add a global feature
new_sys = generic_system(old_sys; energy=10.0u"eV")

# Copy system and add/change cell
new_sys = generic_system(
    old_sys; 
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"],
    periodicity = (true, true, true)
)
```

### Create a Subsystem from an Existing System

The created subsystem does not share any data with the system it was build

```julia
# Subsystem with atoms 1:5
# generic_system(sys, subsys_definition...; kwargs...)
sub_sys = generic_system(sys, 1:5)

# Subsystem with all H and O atoms from sys
sub_sys = generic_system(sys, ChemicalSpecies(:H), ChemicalSpecies(:O))

# Add global feature to subsystem
sub_sys = generic_system(sys, 1:5; label="the first 5 atoms")
```

### Subsystem Views

You can create subsystems that share all data with the host system by calling `system_view`

```julia
# Subsystem with atoms 1:5
# system_view(sys, subsys_definition...)
syb_sys = system_view(sys, 1:5)

# Subsystem with all H and O atoms from sys
sub_sys = system_view(sys, ChemicalSpecies(:H), ChemicalSpecies(:O))
```

Any changes you make to `system_view` structures is made to host system and vise versa.

Note that `system_view` does not see the global features of the host system.

## Changing Systems

AtomsBase defines funtions to modify structures, the following list is supported

- `set_position!(system, i, x)` - all structures
- `set_velocity!(system, i, v)` - all structures that have velocity
- `set_species!(system, i, spc)` - all structures
- `set_cell!(system, cell)` - only for structures with `PeriodicCell` and to another `PeriodicCell` with same dimension. System view structures do not support cell update.
- `set_cell_vectors!(system, bb)` - same as for `set_cell`
- `set_periodicity!(cell, pbc)` - same as for `set_cell`
- `append!(system1, system2)` - if systems have same information fields (e.g. both have velocity), same cell and dont have global features.

**Example**

```julia
sys = generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"],
    [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"];
)

AtomsBase.set_position!(sys, 2, [0.0, 1.0, 0.0]u"Å")
AtomsBase.set_velocity!(sys, 2, [0.0, 0.2, 0.0]u"Å/s")
AtomsBase.set_species!(sys, 1, ChemicalSpecies(:N))
```

## Quality of Life Extensions to AtomsBase

**Methods to change positions**

```julia
center_of_mass(sys)

# translate the whole system by r
translate_system!(sys, r)

# translate a copy of system by r
translate_system(sys, r)

# rotate the system
using Rotations
rot = rand(RotMatrix{3})
rotate_system!(sys, rot)

# rotate a copy of the system
rotate_system(sys, rot)
```

**Add system together or repeat them**

```julia
# make system that has sys2 added to sys1, keep sys1 and sys2 as they are
add_systems(sys1, sys2)

# repeat system along cell vectors
# repeat system 3 times along all cell vectors
repeat(sys, 3)

# repeat 2 times on the first cell vector, 3 times on th esecond cell vector
# and 4 times along the third cell vector
repeat(sys, (2,3,4))
```

**Methods to get information from systems**

```julia
# distance of atoms i and j as a vector
distance_vector(sys, i , j)

# distance of atoms i and j
distance(sys, i, j)

# bond angle of atom i, j and k (j->i vs j->k)
bond_angle(sys, i, j, k)

# dihedral angle of atoms i, j, k and m
dihedral_angle(sys, i, j, k, m)
```

### Fractional Coordinate Methods

```julia
# get inverse cell of the system as matrix
inv_cell(sys)

# get cell_vectors as matrix
cell_matrix(sys)

# fractional coordinates of atom(s) i
fractional_coordinates(sys, i)

# fractional coordinates as matrix for atom(s) i
fractional_coordinates_as_matrix(sys, i)

# wrap atoms inside the cell
wrap_coordinates!(sys)
```