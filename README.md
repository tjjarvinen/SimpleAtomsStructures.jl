# SimpleAtomsStructures

[![Build Status](https://github.com/tjjarvinen/SimpleAtomsStructures.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tjjarvinen/SimpleAtomsStructures.jl/actions/workflows/CI.yml?query=branch%3Amain)


New and improved implementation of AtomsBase interface.

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

# System Structures

Under the hood there are several different system structures that differ on what data they store. But you only need to call `generic_system` to build all of them.
Depending on what information you provide you get different structure

```julia
# Build based on vector of atoms
# generic_system(atoms::AbstractVector{SimpeAtom}; kwargs...)
sys = generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"), 
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]
)

# Same but added key for energy
sys = generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")];
    energy = 10.0u"eV"
)

# Add cell to the system
sys = generic_system([
    SimpleAtom(:H, [0.0, 0.0, 0.0]u"Å"),
    SimpleAtom(:O, [1.0, 0.0, 0.0]u"Å")]; 
    cell_vectors = [[1.0, 0.0, 0.0]u"Å", [0.0, 1.0, 0.0]u"Å", [0.0, 0.0, 1.0]u"Å"], 
    periodicity = (true, true, true)
)

# Create a system from an array of pairs
# generic_system(AbstractVector(<:Pair); kwargs...)
sys = generic_system([:H => [0.0, 0.0, 0.0]u"Å", :O => [1.0, 0.0, 0.0]u"Å"])

# Create a system vectors of atom symbols, positions and velocities
# geric_system(spc, pos, vel)
sys = generic_system(
    [:H, :O],
    [[0.0, 0.0, 0.0]u"Å", [1.0, 0.0, 0.0]u"Å"],
    [[0.1, 0.0, 0.0]u"Å/s", [0.2, 0.0, 0.0]u"Å/s"];
)
```