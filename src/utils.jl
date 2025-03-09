

"""
    center_of_mass(sys)

Return the center of mass of the system `sys`.
"""
function center_of_mass(sys)
    r   = position(sys, :)
    m   = mass(sys, :)
    cms = sum( m.*r ) / sum(m)
    return cms
end

"""
    inv_cell(sys)

Return the inverse of the cell matrix of the AbstractSystem `sys` or cell.
"""
function inv_cell(sys)
    abc = cell_vectors(sys)
    # ignore dimensions with infinities in the cell vectors
    # by replacing them with the unit vector
    abc = map( enumerate(abc) ) do (i, v)
        if all( x -> x.val != Inf, v )
            return v
        else
            tmp = zeros(eltype(v), length(v))
            tmp[i] = one( Unitful.numtype( eltype(v) ) ) * unit( eltype(v) )
            return SVector(tmp...)
        end
    end
    abc_matrix = reduce(hcat, abc)
    abc_inv = inv( ustrip.(abc_matrix) ) / unit(abc_matrix[1,1])
    return abc_inv
end

function inv_cell(::IsolatedCell{D, T}) where{D,T}
    return Matrix{T}(I, D, D)
end

inv_cell(sys::AbstractSystem) = inv_cell(cell(sys))


"""
    cell_matrix(sys)

Return the cell matrix of the `sys` with the cell vectors as columns.

Works when `sys` has `cell_vectors` defined.
"""
function cell_matrix(sys)
    abc = cell_vectors(sys)
    return reduce(hcat, abc)
end

"""
    fractional_coordinates_as_matrix(sys, i)

Return the fractional coordinates of atom(s) `i` in the system `sys` as a matrix.
"""
function fractional_coordinates_as_matrix(sys::AbstractSystem, i)
    abc_inv = inv_cell(sys)
    r = position_as_matrix(sys, i)
    return abc_inv * r
end

fractional_coordinates_as_matrix(sys::AbstractSystem, i::Int) =
    fractional_coordinates(sys, i)

function fractional_coordinates(sys::AbstractSystem, i::Int)
    abc_inv = inv_cell(sys)
    r = position(sys, i)
    return SVector( (abc_inv * r)... )
end

function fractional_coordinates(sys::AbstractSystem{D}, i) where{D}
    tmp = fractional_coordinates_as_matrix(sys, i)
    return reinterpret(reshape, SVector{D, eltype(tmp)}, tmp)
end

function fractional_coordinates_as_matrix(cell::PeriodicCell{D}, coord::AbstractVector{SVector{D,T}}) where{D,T<:Unitful.Length}
    inv_cell_matrix = inv_cell(cell)
    rf = inv_cell_matrix * reinterpret(reshape, T, coord)
    return rf
end

function fractional_coordinates(cell::PeriodicCell{D}, coord::AbstractArray{SVector{D,T}}) where{D,T<:Unitful.Length}
    s = size(coord)
    rf = fractional_coordinates_as_matrix(cell, vec(coord))
    tmp = reinterpret(reshape, SVector{D,eltype(rf)}, rf)
    return Array( reshape(tmp, s) ) # Clear type a bit
end

function fractional_coordinates(cell::PeriodicCell{D}, coord::SVector{D,T}) where{D,T<:Unitful.Length}
    inv_cell_matrix = inv_cell(cell)
    rf = inv_cell_matrix * coord
    return rf    
end


"""
    position_as_matrix(sys, i)

Return the position of atom(s) `i` in the system `sys` as a matrix.
""" 
function position_as_matrix(sys::AbstractSystem, i)
    tmp = position(sys, i)
    ET = (eltype∘eltype)(tmp)
    return reinterpret(reshape, ET, tmp)
end

"""
    wrap_coordinates!(sys)

Wrap the coordinates of the system to the unit cell.
"""
function wrap_coordinates!(sys)
    @argcheck all( periodicity(sys) )
    wrap(x) = x > 0 ? x-ceil(x)+1 : x-floor(x) # 0 <= x < 1
    rf = fractional_coordinates_as_matrix(sys, :)
    rf .= wrap.(rf)
    abc_matrix = cell_matrix(sys)
    tmp = abc_matrix * rf
    T = SVector{n_dimensions(sys), eltype(abc_matrix)}
    new_r = reinterpret(reshape, T, tmp)
    AtomsBase.set_position!(sys, :, new_r)
    return sys
end

wrap_coordinates!(sys::AbstractIsolatedSystem) = sys

function wrap_coordinates!(cell::PeriodicCell{D}, coord::AbstractVector{SVector{D,T}}) where{D, T<:Unitful.Length}
    @argcheck all( periodicity(cell) )
    wrap(x) = x > 0 ? x-ceil(x)+1 : x-floor(x) # 0 <= x < 1
    frc = fractional_coordinates_as_matrix(cell, coord)
    frc .= wrap.(frc)
    abc_matrix = cell_matrix(cell)
    tmp = abc_matrix * frc
    return reinterpret(reshape, SVector{D,T}, tmp)
end

function wrap_coordinates!(cell::PeriodicCell{D}, coord::AbstractArray{SVector{D,T}}) where{D, T<:Unitful.Length}
    s = size(coord)
    rf = wrap_coordinates!(cell, vec(coord))
    return Array( reshape(rf, s) )  # Clear type a bit
end

function wrap_coordinates!(cell::PeriodicCell{D}, coord::SVector{D,T}) where{D,T<:Unitful.Length}
    return wrap_coordinates!(cell, [coord])[1]
end


function translate_system!(sys::AbstractSystem{D}, r::SVector{D, <:Unitful.Length}) where{D}
    for i in 1:length(sys)
        AtomsBase.set_position!(
            sys, i,
            position(sys, i) + r
        )
    end
    return sys
end

"""
    translate_system!(sys, r)

Translate the system `sys` by the vector `r`.
"""
function translate_system!(sys::AbstractSystem{D}, r::AbstractVector{<:Unitful.Length}) where{D}
    @argcheck length(r) == D
    translate_system!(sys, SVector(r...))
    return sys
end

function Base.:+(sys::AbstractSystem{D}, r::SVector{D, <:Unitful.Length}) where{D}
    tmp = deepcopy(sys)
    return translate_system!(tmp, r)
end

function Base.:+(sys::AbstractSystem{D}, r::AbstractVector{<:Unitful.Length}) where{D}
    @argcheck length(r) == D
    return +(sys, SVector(r...))
end

Base.:+(r::AbstractVector{<:Unitful.Length}, sys::AbstractSystem{D}) where{D} = +(sys, r)
Base.:-(sys::AbstractSystem{D}, r::SVector{D, <:Unitful.Length}) where{D} = +(sys, -r)

"""
    rotate_system!(sys, r)

Rotate the system `sys` by the rotation `r`.

Note, this function does work only for isolated systems.
"""
function rotate_system!(sys::AbstractIsolatedSystem, r::Rotation)
    # does not work for all systems in general (e.g. FlexibleSystem)
    pos = position_as_matrix(sys, :)
    pos .= r * pos
    return sys
end

function rotate_system!(sys::GeneralSystem, r::Rotation)
    rotate_system!(sys.base_system, r)
    return sys
end

function Base.:*(sys::AbstractSystem, r::Rotation)
    tmp = deepcopy(sys)
    return rotate_system!(tmp, r)
end

Base.:*(r::Rotation, sys::AbstractSystem) = *(sys, r)


## 

"""
    distance_vector(sys, i, j)

Return the distance vector between atom `i` and atom `j` in the system `sys`.

Note, currently only works for orthorombic cells or isolated systems.
"""
function distance_vector(sys::Union{AbstractIsolatedSystem, AtomsVector}, i::Int, j::Int)
    r1 = position(sys, i)
    r2 = position(sys, j)
    return r2 - r1
end


function distance_vector(sys::AbstractSystem, i::Int, j::Int)
    r1 = fractional_coordinates(sys, i)
    r2 = fractional_coordinates(sys, j)
    return distance_vector(cell(sys), r2 - r1)
end

distance_vector(::IsolatedCell{D}, r::SVector{D}) where{D} = r

function distance_vector(cell::PeriodicCell{D}, r::SVector{D}) where{D}
    wrap(x) = x > 0 ? x-ceil(x)+1 : x-floor(x) # 0 <= x < 1
    function dwrap(x) # distance wrap. x ∈ [-0.5, 0.5]
        if x > 0.5
            return x - 1
        elseif x < -0.5
            return x + 1
        else
            return x
        end
    end
    Δr = map( r, periodicity(cell) ) do x, p
        if p
            return dwrap( wrap(x) )
        else
            return x
        end
    end
    return cell_matrix(cell) * Δr
end


"""
    distance(sys, i, j)

Return the distance between atom `i` and atom `j` in the system `sys`.
"""
function distance(sys::Union{AbstractSystem, AtomsVector}, i::Int, j::Int)
    r = distance_vector(sys, i, j)
    return norm(r)
end

function distance(
    sys1::Union{AbstractIsolatedSystem, AtomsVector},
    sys2::Union{AbstractIsolatedSystem, AtomsVector}
)
    return [ norm( position(sys2, j) - position(sys1,i) ) for i in 1:length(sys1), j in 1:length(sys2) ]
end

"""
    angle(sys, i, j, k)

Calculate the angle between atoms `i`, `j`, and `k` in the system `sys`.

The angle is defined as the angle between the vectors `r_ij` (from atom `j` to atom `i`)
and `r_jk` (from atom `j` to atom `k`).

See also `angled`.
"""
function Base.angle(sys, i::Int, j::Int, k::Int)
    r1 = distance_vector(sys, j, i)
    r2 = distance_vector(sys, j, k)
    return acos(dot(r1,r2)/sqrt(dot(r1,r1)*dot(r2,r2)))
end

"""
    angled(sys, i, j, k)

Calculate the angle between atoms `i`, `j`, and `k` in the system `sys`.

The angle is defined as the angle between the vectors `r_ij` (from atom `j` to atom `i`)
and `r_jk` (from atom `j` to `k`). The result is in degrees.

See also `angle`.
"""
function angled(sys, i::Int, j::Int, k::Int)
    r1 = distance_vector(sys, j, i)
    r2 = distance_vector(sys, j, k)
    return acosd(dot(r1,r2)/sqrt(dot(r1,r1)*dot(r2,r2)))
end

"""
    dihedral_angle(sys, i, j, k, m)

Calculate the dihedral angle between atoms `i`, `j`, `k`, and `m` in the system `sys`.

The dihedral angle is defined as the angle between the planes defined by the atoms `i`, `j`, `k`
and `j`, `k`, `m`.

See also `dihedral_angled`.
"""
function dihedral_angle(sys, i::Int, j::Int, k::Int, m::Int)
    r1 = distance_vector(sys, i,j)
    r2 = distance_vector(sys, j,k)
    r3 = distance_vector(sys, k,m)
    t1 = cross(cross(r1,r2), cross(r2,r3))
    t2 = dot(cross(r1,r2), cross(r2,r3))
    return atan( dot(t1, r2./norm(r2)), t2 )
end


"""
    dihedral_angled(sys, i, j, k, m)

Calculate the dihedral angle between atoms `i`, `j`, `k`, and `m` in the system `sys`.

The dihedral angle is defined as the angle between the planes defined by the atoms `i`, `j`, `k`
and `j`, `k`, `m`. The result is in degrees.

See also `dihedral_angle`.
"""
function dihedral_angled(sys, i::Int, j::Int, k::Int, m::Int)
    r1 = distance_vector(sys, i,j)
    r2 = distance_vector(sys, j,k)
    r3 = distance_vector(sys, k,m)
    t1 = cross(cross(r1,r2), cross(r2,r3))
    t2 = dot(cross(r1,r2), cross(r2,r3))
    return atand( dot(t1, r2./norm(r2)), t2 )
end

##

# make this with generating function for more dimensions
function Base.repeat(sys::CellSystem{3}, n::NTuple{3,<:Integer}) # where{D}
    abc = cell_vectors(sys)
    abc_n = n .* abc
    cell = PeriodicCell(abc_n, periodicity(sys))
    tsys = CellSystem(sys.base_system, cell)
    nsys = deepcopy(tsys)
    for i in 0:n[1]-1, j in 0:n[2]-1, k in 0:n[3]-1
        if ! ( 0 == i == j == k )
            r = sum( [i,j,k] .* abc )
            tmp = tsys + r
            append!(nsys, tmp)
        end
    end
    return nsys
end

Base.repeat(sys::CellSystem{3}, n::Integer) = Base.repeat(sys, (n,n,n))

Base.:*(sys::CellSystem{3}, n::Integer) = Base.repeat(sys, n)
Base.:*(sys::CellSystem{3}, n::NTuple{3,<:Integer}) = Base.repeat(sys, n)
Base.:*(n::Integer, sys::CellSystem{3}) = Base.repeat(sys, n)
Base.:*(n::NTuple{3,<:Integer}, sys::CellSystem{3}) = Base.repeat(sys, n)