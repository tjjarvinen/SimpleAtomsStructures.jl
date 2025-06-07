using AtomsBase
using AtomsBaseTesting
using SimpleAtomsStructures
using Rotations
using Unitful
using Test

@testset "SimpleAtomsStructures.jl" begin
    # Write your tests here.
    ref = make_test_system()
    @testset "SimpleSystem" begin
        sys = SimpleSystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy] 
    end
    @testset "SimpleVelocitySystem" begin
        sys = SimpleVelocitySystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy]
    end
    @testset "AtomicPropertySystem" begin
        sys = AtomicPropertySystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test all( mass(sys, :) .== mass( ref.system, :) )
        @test isa(cell(sys), IsolatedCell)
        @test all( sys[:periodicity] .== (false, false, false) )
        @test_throws KeyError sys[:dummy]
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
    end
    @testset "CellSystemSystem" begin
        sys = CellSystem(ref.system)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test all( mass(sys, :) .== mass( ref.system, :) )
        @test isa(cell(sys), PeriodicCell)
        @test all( cell_vectors(sys) .≈ ref.cell_vectors )
        @test all( sys[:periodicity] .== ref.periodicity )
        @test all( sys[:cell_vectors] .≈ ref.cell_vectors )
        @test_throws KeyError sys[:dummy]
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
    end
    @testset "GeneralSystem" begin
        sys = GenericSystem(ref.system)
        @test isa(sys, SimpleAtomsStructures.GeneralSystem)
        @test all( position(sys, :) .≈ position( ref.system, :) )
        @test all( velocity(sys, :) .≈ velocity( ref.system, :) )
        @test all( species(sys, :) .== species( ref.system, :) )
        @test all( mass(sys, :) .== mass( ref.system, :) )
        @test isa(cell(sys), PeriodicCell)
        @test all( cell_vectors(sys) .≈ ref.cell_vectors )
        @test all( x-> hasatomkey(sys, x), keys(ref.atprop) )
        @test all( x-> haskey(ref.atprop, x), atomkeys(sys) )
        for (k,v) in pairs(ref.atprop)
            sys_val = [ x[k] for x in sys ]
            @test all( sys_val .== v)
        end
        @test all( x-> haskey(sys, x), keys(ref.sysprop) )
        @test all( x-> haskey(ref.sysprop, x), keys(sys) )
        @test all( [all(sys[k] .== v) for (k,v) in pairs(ref.sysprop)] )
        @test_throws KeyError sys[:dummy]
    end
    @testset "Utils" begin
        sys = SimpleSystem(ref.system)
        # rotation tests
        q = rand(QuatRotation)
        sys2 = rotate_system(sys, q)
        @test all( i-> position(sys2, i) ≈ q * position(sys, i), 1:length(sys) )
        @test angle(sys, 1, 2, 3) ≈ angle(sys2, 1, 2, 3)
        @test dihedral_angle(sys, 1,2,3,4) ≈ dihedral_angle(sys2, 1,2,3,4)
        @test distance_vector(sys2, 1, 2) ≈ q * distance_vector(sys, 1, 2)

        # translation tests
        cms = center_of_mass(sys)
        sys2 = translate_system(sys, -cms)
        @test all( i-> position(sys2, i) ≈ position(sys, i) - cms, 1:length(sys) )
        @test angle(sys, 1, 2, 3) ≈ angle(sys2, 1, 2, 3)
        @test dihedral_angle(sys, 1,2,3,4) ≈ dihedral_angle(sys2, 1,2,3,4)
        @test distance_vector(sys2, 1, 2) ≈ distance_vector(sys, 1, 2)

        
        sys3 = add_systems(sys, sys2)
        @test isa(sys3, typeof(sys))
        @test length(sys3) == length(sys) + length(sys2)
        
        # now with cell
        sys = GenericSystem(ref.system)

        fp = fractional_coordinates(sys, :)
        @test length(fp) == length(sys)
        fpm = fractional_coordinates_as_matrix(sys, :)
        @test size(fpm) == (3, length(sys))

        clm = cell_matrix(sys)
        clv = cell_vectors(sys)
        @test size(clm) == (3,3)
        @test all( x -> all(x[1] .≈ x[2]), zip(clv, eachcol(clm)) )
        icell = inv_cell(sys)
        tmp = icell * clm

        # Repeat system
        csys = CellSystem(ref.system)
        sys123 = repeat(csys, (1, 2, 3))
        @test length(sys123) == 6 * length(csys)
        @test all( species( sys123, :) .=== repeat( species(csys,:), 6) )
        c1 = cell_vectors(csys)
        c2 = cell_vectors(sys123)
        @test c1[1] ≈   c1[1]
        @test c2[2] ≈ 2*c1[2]
        @test c2[3] ≈ 3*c1[3]      
    end
    @testset "SimpleAtom" begin
        sys = GenericSystem(ref.system)
        va = sys[:]
        @test all( k -> k in atomkeys(sys), atomkeys(va) )
        @test all( k -> k in atomkeys(va), atomkeys(sys) )
        @test all( mass(va, :) .≈ mass(sys, :) )
        @test all( species(va, :) .== species(sys, :) )
        @test all( position(va, :) .≈ position(sys, :) )
        @test all( velocity(va, :) .≈ velocity(sys, :) )
    end
    @testset "Trajectory" begin
        sys = GenericSystem(ref.system)
        sys2 = deepcopy( sys )
        translate_system!(sys2, [1., 2., 3.]u"Å")
        trj = VelocityTrajectory([sys, sys2])
        @test length(trj) == 2
        @test all( species(trj, :) .=== species(sys, :) )
        @test all( position(trj, :, 1) .≈ position(sys, :) )
        @test all( position(trj, :, 2) .≈ position(sys2, :) )
        @test all( velocity(trj, :, 1) .≈ velocity(sys, :) )
        @test all( velocity(trj, :, 2) .≈ velocity(sys2, :) )
        @test cell(trj, 1) == cell(sys)
        @test cell(trj, 2) == cell(sys2) 
        @test all( distance(trj, 1, 2, 1:2) .≈ distance(trj, 1, 2, :) )
        @test all( angle(trj, 1, 2, 3, 1:2) .≈ angle(trj, 1, 2, 3, :) )
        @test all( dihedral_angle(trj, 1, 2, 3, 4, 1:2) .≈ dihedral_angle(trj, 1, 2, 3, 4, :) )

        @testset "SystemView" begin
            sv = SystemView(trj, 1)
            @test all( species(sv, :) .=== species(sys, :) )
            @test all( position(sv, :) .≈ position(sys, :) )
            @test all( velocity(sv, :) .≈ velocity(sys, :) )
            @test cell(sv) == cell(sys)
            @test all( periodicity(sv) .== periodicity(sys) )
            @test all( cell_vectors(sv) .≈ cell_vectors(sys) )
            sv2 = SystemView(trj, 2)
            @test all( species(sv2, :) .=== species(sys2, :) )
            @test all( position(sv2, :) .≈ position(sys2, :) )
            @test all( velocity(sv2, :) .≈ velocity(sys2, :) )
            @test cell(sv2) == cell(sys2)

        end

    end
end
