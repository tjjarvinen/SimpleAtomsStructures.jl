

function center_of_mass(sys)
    r   = position(sys, :)
    m   = mass(sys, :)
    cms = sum( m.*r ) / sum(m)
    return cms
end

function inv_cell(sys)
    abc = cell_vectors(sys)
    abc_matrix = reduce(hcat, abc)
    abc_inv = inv( ustrip.(abc_matrix) ) / unit(abc_matrix[1,1])
    return abc_inv
end

function cell_matrix(sys)
    abc = cell_vectors(sys)
    return reduce(hcat, abc)
end

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
    return abc_inv * r
end

function fractional_coordinates(sys::AbstractSystem{D}, i) where{D}
    tmp = fractional_coordinates_as_matrix(sys, i)
    return reinterpret(reshape, SVector{D, eltype(tmp)}, tmp)
end


function position_as_matrix(sys::AbstractSystem, i)
    tmp = position(sys, i)
    ET = (eltype∘eltype)(tmp)
    return reinterpret(reshape, ET, tmp)
end


function wrap_coordinates!(sys)
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


function translate_system!(sys::AbstractSystem{D}, r::SVector{D, <:Unitful.Length}) where{D}
    for i in 1:length(sys)
        AtomsBase.set_position!(
            sys, i,
            position(sys, i) + r
        )
    end
    return sys
end

function Base.:+(sys::AbstractSystem{D}, r::SVector{D, <:Unitful.Length}) where{D}
    tmp = deepcopy(sys)
    return translate_system!(tmp, r)
end

Base.:+(r::SVector{D, <:Unitful.Length}, sys::AbstractSystem{D}) where{D} = +(sys, r)
Base.:-(sys::AbstractSystem{D}, r::SVector{D, <:Unitful.Length}) where{D} = +(sys, -r)


function rotate_system!(sys::AbstractIsolatedSystem, r::Rotation)
    pos = position_as_matrix(sys, :)
    pos .*= r
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

function distance_vector(sys::Union{AbstractIsolatedSystem, AtomsVector}, i::Int, j::Int)
    r1 = position(sys, i)
    r2 = position(sys, j)
    return r2 - r1
end


function distance_vector(sys::CellSystem, i::Int, j::Int)
    # only works for orthorombic cells
    @argcheck periodicity(sys) == (true, true, true)
    wrap(x) = x > 0 ? x-ceil(x)+1 : x-floor(x) # 0 <= x < 1
    function dwrap(x) # distance wrap
        if x > 0.5
            return x - 1
        elseif x < -0.5
            return x + 1
        else
            return x
        end
    end
    r1 = wrap.( fractional_coordinates(sys, i) )
    r2 = wrap.( fractional_coordinates(sys, j) )
    tmp = dwrap.( r2 - r1 )
    cell = cell_matrix(sys)
    return cell * tmp
end

function distance_vector(sys::GeneralSystem, i::Int, j::Int)
    return distance_vector(sys.base_system, i ,j)
end


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

function angle(sys, i::Int, j::Int, k::Int)
    r1 = distance_vector(sys, j, i)
    r2 = distance_vector(sys, j, k)
    return acos(dot(r1,r2)/sqrt(dot(r1,r1)*dot(r2,r2)))
end

function angled(sys, i::Int, j::Int, k::Int)
    r1 = distance_vector(sys, j, i)
    r2 = distance_vector(sys, j, k)
    return acosd(dot(r1,r2)/sqrt(dot(r1,r1)*dot(r2,r2)))
end


function dihedral_angle(sys, i::Int, j::Int, k::Int, m::Int)
    r1 = distance_vector(sys, i,j)
    r2 = distance_vector(sys, j,k)
    r3 = distance_vector(sys, k,m)
    t1 = cross(cross(r1,r2), cross(r2,r3))
    t2 = dot(cross(r1,r2), cross(r2,r3))
    return atan( dot(t1, r2./norm(r2)), t2 )
end

function dihedral_angled(sys, i::Int, j::Int, k::Int, m::Int)
    r1 = distance_vector(sys, i,j)
    r2 = distance_vector(sys, j,k)
    r3 = distance_vector(sys, k,m)
    t1 = cross(cross(r1,r2), cross(r2,r3))
    t2 = dot(cross(r1,r2), cross(r2,r3))
    return atand( dot(t1, r2./norm(r2)), t2 )
end

##


function Base.repeat(sys::CellSystem{3}, n::NTuple{3,Int}) # where{D}
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

Base.repeat(sys::CellSystem{3}, n::Int) = Base.repeat(sys, (n,n,n))

Base.:*(sys::CellSystem{3}, n::Int) = Base.repeat(sys, n)
Base.:*(sys::CellSystem{3}, n::NTuple{3,Int}) = Base.repeat(sys, n)
Base.:*(n::Int, sys::CellSystem{3}) = Base.repeat(sys, n)
Base.:*(n::NTuple{3,Int}, sys::CellSystem{3}) = Base.repeat(sys, n)