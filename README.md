# SimpleAtomsStructures

[![Build Status](https://github.com/tjjarvinen/SimpleAtomsStructures.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tjjarvinen/SimpleAtomsStructures.jl/actions/workflows/CI.yml?query=branch%3Amain)


Design idea here is that there is a simple structure with only species and position information. This structure is then wrapped to other structures that add more information 

You will have either 
1. `SimpleSystem` -> `AtomicPropertySystem` -> `CellSystem` -> `GeneralSystem`

or 

2. `SimpleVelocitySystem` -> `AtomicPropertySystem` -> `CellSystem` -> `GeneralSystem`

depending on what information you need to store.

All systems are build with `GenericSysten` call

```julia
GenericSystem(AbstractSystem, <optional subsystem definition>; global_features...)
GenericSystem(species, positions; global_features...)
GenericSystem(species, positions, velocity; global_features...)
GenericSystem(AtomsVector; global_features...)
```

## Different Systems Structures

Initilize with

```julia
using AtomsBase
using AtomsBaseTesting
using SimpleAtomsStructures


ref = make_test_system()
```

After that you can create systems by

```julia

# Only species and positions
julia> SimpleSystem(ref.system)
SimpleSystem(CHeH₂N, periodicity = FFF):
    SimpleAtom(H,  [-1.19999, -0.521534,  1.17446]u"Å")
    SimpleAtom(H,  [0.532801, -0.653995, 0.145677]u"Å")
    SimpleAtom(C,  [ 1.08169, -0.627532,  1.37471]u"Å")
    SimpleAtom(N,  [ 2.85773, -0.37975, 0.572554]u"Å")
    SimpleAtom(He, [-2.00971, -0.875141,  1.14164]u"Å")

# Add velocities
julia> SimpleVelocitySystem(ref.system)
SimpleVelocitySystem(CHeH₂N, periodicity = FFF):
    SimpleAtom(H,  [-1.19999, -0.521534,  1.17446]u"Å", [ -306058,  -660815, 1.51741e+06]u"m s^-1")
    SimpleAtom(H,  [0.532801, -0.653995, 0.145677]u"Å", [  498585,   254433,  -798499]u"m s^-1")
    SimpleAtom(C,  [ 1.08169, -0.627532,  1.37471]u"Å", [-17145.5,  -835281,  -645397]u"m s^-1")
    SimpleAtom(N,  [ 2.85773, -0.37975, 0.572554]u"Å", [1.08308e+06,  -713829,  -420172]u"m s^-1")
    SimpleAtom(He, [-2.00971, -0.875141,  1.14164]u"Å", [-2.34931e+06,   197953,   597271]u"m s^-1")

# Add general atom features
julia> sys = AtomicPropertySystem(ref.system);
julia> atomkeys(sys)
(:position, :velocity, :species, :mass, :magnetic_moment, :charge, :vdw_radius, :covalent_radius)

# All but global features
julia> CellSystem(ref.system)
CellSystem(CHeH₂N, periodicity = TTF):
    cell_vectors      : [ 1.50304 0.850344 0.717239;
                          0.36113  1.00814 0.814712;
                          0.06828 0.381122  2.12908]u"Å"

    SimpleAtom(H,  [-1.19999, -0.521534,  1.17446]u"Å", [ -306058,  -660815, 1.51741e+06]u"m s^-1")
    SimpleAtom(H,  [0.532801, -0.653995, 0.145677]u"Å", [  498585,   254433,  -798499]u"m s^-1")
    SimpleAtom(C,  [ 1.08169, -0.627532,  1.37471]u"Å", [-17145.5,  -835281,  -645397]u"m s^-1")
    SimpleAtom(N,  [ 2.85773, -0.37975, 0.572554]u"Å", [1.08308e+06,  -713829,  -420172]u"m s^-1")
    SimpleAtom(He, [-2.00971, -0.875141,  1.14164]u"Å", [-2.34931e+06,   197953,   597271]u"m s^-1")

# Everything supported
julia> GenericSystem(ref.system)  # is equal to ref.system
GeneralSystem(CHeH₂N, periodicity = TTF):
    cell_vectors      : [ 1.50304 0.850344 0.717239;
                          0.36113  1.00814 0.814712;
                          0.06828 0.381122  2.12908]u"Å"
    multiplicity      : 2
    charge            : -1 e
    extra_data        : 42

    SimpleAtom(H,  [-1.19999, -0.521534,  1.17446]u"Å", [ -306058,  -660815, 1.51741e+06]u"m s^-1")
    SimpleAtom(H,  [0.532801, -0.653995, 0.145677]u"Å", [  498585,   254433,  -798499]u"m s^-1")
    SimpleAtom(C,  [ 1.08169, -0.627532,  1.37471]u"Å", [-17145.5,  -835281,  -645397]u"m s^-1")
    SimpleAtom(N,  [ 2.85773, -0.37975, 0.572554]u"Å", [1.08308e+06,  -713829,  -420172]u"m s^-1")
    SimpleAtom(He, [-2.00971, -0.875141,  1.14164]u"Å", [-2.34931e+06,   197953,   597271]u"m s^-1")

# same as useing isolated_system from AtomsBase
julia> hydrogen = GenericSystem([ 
           :H => [0.1, 0, 0.]u"Å",
           :H => [0, 0, 1.]u"Å",
           :H => [4., 0, 0.]u"Å",
           :H => [4., 1., 0.]u"Å"
       ])
SimpleSystem(H₄, periodicity = FFF):
    SimpleAtom(H,  [     0.1,        0,        0]u"Å")
    SimpleAtom(H,  [       0,        0,        1]u"Å")
    SimpleAtom(H,  [       4,        0,        0]u"Å")
    SimpleAtom(H,  [       4,        1,        0]u"Å")
```

## Changing structures

```julia
# Add global features
julia> GenericSystem(ref.system; my_property="hello")
GeneralSystem(CHeH₂N, periodicity = TTF):
    cell_vectors      : [ 1.50304 0.850344 0.717239;
                          0.36113  1.00814 0.814712;
                          0.06828 0.381122  2.12908]u"Å"
    multiplicity      : 2
    charge            : -1 e
    my_property       : hello
    extra_data        : 42

    SimpleAtom(H,  [-1.19999, -0.521534,  1.17446]u"Å", [ -306058,  -660815, 1.51741e+06]u"m s^-1")
    SimpleAtom(H,  [0.532801, -0.653995, 0.145677]u"Å", [  498585,   254433,  -798499]u"m s^-1")
    SimpleAtom(C,  [ 1.08169, -0.627532,  1.37471]u"Å", [-17145.5,  -835281,  -645397]u"m s^-1")
    SimpleAtom(N,  [ 2.85773, -0.37975, 0.572554]u"Å", [1.08308e+06,  -713829,  -420172]u"m s^-1")
    SimpleAtom(He, [-2.00971, -0.875141,  1.14164]u"Å", [-2.34931e+06,   197953,   597271]u"m s^-1")

# Take a subsystem (and discard global features)
julia> GenericSystem(ref.system, 1:3; subs_name="my subsystem")
GeneralSystem(CH₂, periodicity = TTF):
    cell_vectors      : [ 1.50304 0.850344 0.717239;
                          0.36113  1.00814 0.814712;
                          0.06828 0.381122  2.12908]u"Å"
    subs_name         : my subsystem

    SimpleAtom(H,  [-1.19999, -0.521534,  1.17446]u"Å", [ -306058,  -660815, 1.51741e+06]u"m s^-1")
    SimpleAtom(H,  [0.532801, -0.653995, 0.145677]u"Å", [  498585,   254433,  -798499]u"m s^-1")
    SimpleAtom(C,  [ 1.08169, -0.627532,  1.37471]u"Å", [-17145.5,  -835281,  -645397]u"m s^-1")

# Subsystem with only Hydrogen and Helium
julia> GenericSystem(ref.system, ChemicalSpecies(:H), ChemicalSpecies(:He); id=1)
GeneralSystem(HeH₂, periodicity = TTF):
    cell_vectors      : [ 1.50304 0.850344 0.717239;
                          0.36113  1.00814 0.814712;
                          0.06828 0.381122  2.12908]u"Å"
    id                : 1

    SimpleAtom(H,  [-1.19999, -0.521534,  1.17446]u"Å", [ -306058,  -660815, 1.51741e+06]u"m s^-1")
    SimpleAtom(H,  [0.532801, -0.653995, 0.145677]u"Å", [  498585,   254433,  -798499]u"m s^-1")
    SimpleAtom(He, [-2.00971, -0.875141,  1.14164]u"Å", [-2.34931e+06,   197953,   597271]u"m s^-1")
```

## Atom structure

There is a `SimpleAtom` structure that represents atoms

```julia
julia> SimpleAtom(:H, zeros(3)*u"Å")
SimpleAtom(H, Z = 1, m = 1.008 u):
    position          : [0,0,0]u"Å"
    species           : H

# Add extra properties
julia> atoms = [
          SimpleAtom(:H, zeros(3)*u"Å"; my_property=1),
          SimpleAtom(:H, ones(3)*u"Å"; my_property=2)
       ]
2-element Vector{SimpleAtom{3, @NamedTuple{species::ChemicalSpecies, position::StaticArraysCore.SVector{3, Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}, my_property::Int64}, Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}}}:
 SimpleAtom(H,  [       0,        0,        0]u"Å")
 SimpleAtom(H,  [       1,        1,        1]u"Å")

# Vector of SimpleAtoms has core AtomsBase interface implemented
julia> position(atoms, 2)
3-element StaticArraysCore.SVector{3, Quantity{Float64, 𝐋, Unitful.FreeUnits{(Å,), 𝐋, nothing}}} with indices SOneTo(3):
 1.0 Å
 1.0 Å
 1.0 Å

# Create full system from vector of atoms
julia> GenericSystem(atoms)
AtomicPropertySystem(H₂, periodicity = FFF):
    SimpleAtom(H,  [       0,        0,        0]u"Å")
    SimpleAtom(H,  [       1,        1,        1]u"Å")



```