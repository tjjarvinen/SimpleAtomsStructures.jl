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

function _load_cell(frame::Chemfiles.Frame)
    ccell = Chemfiles.UnitCell(frame)
    cell_shape = Chemfiles.shape(ccell)
    if cell_shape == Chemfiles.Infinite
        return IsolatedCell(3)
    else
        cell_vectors = collect(eachrow(Chemfiles.matrix(ccell)))u"Å"
        return PeriodicCell(; cell_vectors, periodicity=(true, true, true))
    end 
end

function SimpleAtomsStructures.CellSystem(frame::Chemfiles.Frame)
    sys = SimpleAtomsStructures.SimpleVelocitySystem(frame)
    cell = _load_cell(frame)
    return SimpleAtomsStructures.CellSystem(sys, cell)
end

function SimpleAtomsStructures.SimpleTrajectory(traj::Chemfiles.Trajectory)
    # SimpleTrajectory has constant cell, so we can use the first frame to initialize the system
    first_frame = Chemfiles.read_step(traj, 0)
    sys = SimpleAtomsStructures.SimpleSystem(first_frame)
    ntraj = SimpleAtomsStructures.SimpleTrajectory(sys)
    for frame in traj
        pos = positions(frame) * u"Å"
        r = reinterpret(reshape, SVector{3, Float64}, pos) * u"Å"
        append!(ntraj, r)
    end
    return ntraj 
end


function SimpleAtomsStructures.SimpleVelocityTrajectory(traj::Chemfiles.Trajectory)
    # SimpleVelocityTrajectory has constant cell, so we can use the first frame to initialize the system
    first_frame = Chemfiles.read_step(traj, 0)
    if has_velocities(first_frame)
        sys = SimpleAtomsStructures.SimpleVelocitySystem(first_frame)
        ntraj = SimpleAtomsStructures.SimpleVelocityTrajectory(sys)
        for frame in traj
            pos = positions(frame) * u"Å"
            vel = velocities(frame) * u"Å/ps"
            r = reinterpret(reshape, SVector{3, Float64}, pos) * u"Å"
            v = reinterpret(reshape, SVector{3, Float64}, vel) * u"Å/ps"
            append!(ntraj, r, v)
        end
        return ntraj
    end
    return SimpleAtomsStructures.SimpleTrajectory(traj)
end

function SimpleAtomsStructures.ConstantVolumeTrajectory(traj::Chemfiles.Trajectory)
    # ConstantVolumeTrajectory has constant cell, so we can use the first frame to initialize the system
    first_frame = Chemfiles.read_step(traj, 0)
    tcell = _load_cell(first_frame)
    btraj = SimpleAtomsStructures.SimpleVelocityTrajectory(traj)
    return SimpleAtomsStructures.ConstantVolumeTrajectory(btraj, tcell)
end

function SimpleAtomsStructures.SimpleTrajectory(fname::AbstractString)
    traj = Chemfiles.Trajectory(fname)
    return SimpleAtomsStructures.SimpleTrajectory(traj)
end

function SimpleAtomsStructures.SimpleVelocityTrajectory(fname::AbstractString)
    traj = Chemfiles.Trajectory(fname)
    return SimpleAtomsStructures.SimpleVelocityTrajectory(traj)    
end

function SimpleAtomsStructures.ConstantVolumeTrajectory(fname::AbstractString; species_from=nothing)
    ttraj = Chemfiles.Trajectory(fname)
    traj = SimpleAtomsStructures.ConstantVolumeTrajectory(ttraj)
    if isnothing(species_from)
        return traj
    end
    tmp_traj = Chemfiles.Trajectory(species_from)
    frame = Chemfiles.read(tmp_traj)
    tmp = SimpleAtomsStructures.SimpleSystem(frame)
    sys = traj[1]
    if length(tmp) == length(sys)
        AtomsBase.set_species!(sys, :, species(tmp, :))
        return traj
    end
    error("Species vector length does not match the number of atoms in the trajectory")
end

end