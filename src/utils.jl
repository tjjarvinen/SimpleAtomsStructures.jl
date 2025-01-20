

function center_of_mass(sys)
    r   = position(sys, :)
    m   = mass(sys, :)
    cms = sum( m.*r ) / sum(m)
    return cms
end

function fractional_coordinates(sys, i)
    abc = cell_vectors(sys)
    abc_matrix = reduce(hcat, abc)
    abc_inv = inv( ustrip.(abc_matrix) )
    r = position(sys, i)
    T = Unitful.numtype( (eltype∘eltype)(r) )
    r_matrix = reinterpret(reshape, T, r )
    return abc_inv * r_matrix
end


function wrap_coordinates!(sys)
    wrap(x) = x > 0 ? x-ceil(x)+1 : x-floor(x) # 0 <= x < 1
    rf = fractional_coordinates(sys, :)
    rf .= wrap.(rf)
    abc = cell_vectors(sys)
    abc_matrix = reduce(hcat, abc)
    tmp = abc_matrix * rf
    T = SVector{n_dimensions(sys), eltype(abc_matrix)}
    new_r = reinterpret(reshape, T, tmp)
    AtomsBase.set_position!(sys, :, new_r)
    return sys
end


function translate!(sys::AbstractSystem{D}, r::SVector{D, <:Unitful.Length}) where{D}
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
    return translate!(tmp, r)
end

Base.:+(r::SVector{D, <:Unitful.Length}, sys::AbstractSystem{D}) where{D} = +(sys, r)
Base.:-(r::SVector{D, <:Unitful.Length}, sys::AbstractSystem{D}) where{D} = +(sys, -r)
Base.:-(sys::AbstractSystem{D}, r::SVector{D, <:Unitful.Length}) where{D} = +(sys, -r)