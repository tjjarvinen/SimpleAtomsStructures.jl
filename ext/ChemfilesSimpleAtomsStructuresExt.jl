module ChemfilesSimpleAtomsStructuresExt

using AtomsBase
using Chemfiles
using SimpleAtomsStructures
using Unitful
using StaticArrays


function SimpleAtomsStructures.SimpleSystem(frame::Chemfiles.Frame)
    pos = positions(frame)
    r = reinterpret(reshape, SVector{3, Float64}, pos) * u"Å"
    spc = map( frame ) do a
        ChemicalSpecies( atomic_number(a) )
    end
    SimpleAtomsStructures.SimpleSystem(spc, r)
end

function SimpleAtomsStructures.SimpleVelocitySystem(frame::Chemfiles.Frame)
    pos = positions(frame)
    r = reinterpret(reshape, SVector{3, Float64}, pos) * u"Å"
    spc = map( frame ) do a
        ChemicalSpecies( atomic_number(a) )
    end
    if has_velocities(frame)
        vel = velocities(frame)
        v = reinterpret(reshape, SVector{3, Float64}, vel) * u"Å/ps"
        return SimpleAtomsStructures.SimpleVelocitySystem(spc, r, v)
    end
    return SimpleAtomsStructures.SimpleSystem(spc, r)
end


function SimpleAtomsStructures.CellSystem(frame::Chemfiles.Frame)
    sys = SimpleAtomsStructures.SimpleVelocitySystem(frame)

    cell_shape = Chemfiles.shape(Chemfiles.UnitCell(frame))
    if cell_shape == Chemfiles.Infinite
        cell = IsolatedCell(3)
    else
        @assert cell_shape in (Chemfiles.Triclinic, Chemfiles.Orthorhombic)
        cell_vectors = collect(eachrow(Chemfiles.matrix(Chemfiles.UnitCell(frame))))u"Å"
        cell = PeriodicCell(; cell_vectors, periodicity=(true, true, true))
    end

    return SimpleAtomsStructures.CellSystem(sys, cell)
end

function SimpleAtomsStructures.SimpleTrajectory(traj::Chemfiles.Trajectory)
    # SimpleTrajectory has constant cell, so we can use the first frame to initialize the system
    first_frame = Chemfiles.read_step(traj, 0)
    sys = CellSystem(first_frame)
    ntraj = SimpleAtomsStructures.SimpleTrajectory(sys)
    for frame in traj
        pos = positions(frame) * u"Å"
        push!(ntraj, pos)
    end
    return ntraj 
end


function SimpleAtomsStructures.SimpleVelocityTrajectory(traj::Chemfiles.Trajectory)
    # SimpleVelocityTrajectory has constant cell, so we can use the first frame to initialize the system
    first_frame = Chemfiles.read_step(traj, 0)
    if !has_velocities(first_frame)
        return SimpleAtomsStructures.SimpleTrajectory(traj)
    end
    sys = CellSystem(first_frame)
    ntraj = SimpleAtomsStructures.SimpleVelocityTrajectory(sys)
    for frame in traj
        pos = positions(frame) * u"Å"
        vel = velocities(frame) * u"Å/ps"
        push!(ntraj, pos, vel)
    end
    return ntraj
end

function SimpleAtomsStructures.SimpleTrajectory(fname::AbstractString)
    traj = Chemfiles.Trajectory(fname)
    return SimpleAtomsStructures.SimpleTrajectory(traj)
end

function SimpleAtomsStructures.SimpleVelocityTrajectory(fname::AbstractString)
    traj = Chemfiles.Trajectory(fname)
    return SimpleAtomsStructures.SimpleVelocityTrajectory(traj)    
end

end