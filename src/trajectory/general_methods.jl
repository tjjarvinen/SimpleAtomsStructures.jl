

function SimpleAtomsStructures.distance(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, frame::Int=1)
    return SimpleAtomsStructures.distance(traj[frame], i, j)
end

function SimpleAtomsStructures.distance(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, ::Colon)
    map(traj) do sys
        SimpleAtomsStructures.distance(sys, i, j)
    end
end

function SimpleAtomsStructures.distance(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, frames)
    map( frames ) do frame
        SimpleAtomsStructures.distance(traj[frame], i, j)
    end
end

function SimpleAtomsStructures.dihedral_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, l::Int, frame::Int)
    return SimpleAtomsStructures.dihedral_angle(traj[frame], i, j, k, l)
end

function SimpleAtomsStructures.dihedral_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, l::Int, ::Colon)
    map(traj) do sys
        SimpleAtomsStructures.dihedral_angle(sys, i, j, k, l)
    end
end

function SimpleAtomsStructures.dihedral_angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, l::Int, frames)
    map( frames ) do frame
        SimpleAtomsStructures.dihedral_angle(traj[frame], i, j, k, l)
    end
end


function Base.angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, frame::Int)
    return SimpleAtomsStructures.angle(traj[frame], i, j, k)
end

function Base.angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, ::Colon)
    map(traj) do sys
        SimpleAtomsStructures.angle(sys, i, j, k)
    end     
end

function Base.angle(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, k::Int, frames)
    map( frames ) do frame
        SimpleAtomsStructures.angle(traj[frame], i, j, k)
    end       
end

function SimpleAtomsStructures.distance_vector(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, frame::Int)
    return SimpleAtomsStructures.distance_vector(traj[frame], i, j)
end

function SimpleAtomsStructures.distance_vector(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, ::Colon)
    map(traj) do sys
        SimpleAtomsStructures.distance_vector(sys, i, j)
    end
end

function SimpleAtomsStructures.distance_vector(traj::AbstractVector{<:AbstractSystem}, i::Int, j::Int, frames)
    map( frames ) do frame
        SimpleAtomsStructures.distance_vector(traj[frame], i, j)
    end
end