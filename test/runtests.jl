using Test
using CleanroomSolver
using CleanroomSolver.Common
using CleanroomSolver.Grid
using CleanroomSolver.Fields
using CleanroomSolver.Convection
using CleanroomSolver.PressureSolver
using CleanroomSolver.FractionalStep
using CleanroomSolver.BoundaryConditions
using CleanroomSolver.InputReader

@testset "CleanroomSolver Tests" begin

    @testset "1. Grid Generation" begin
        # Setup config
        nx, ny, nz = 10, 10, 10
        dim = DimensionParams(1.0, 1.0, 0.01, 100.0, 1.0)
        conf = GridConfig(nx, ny, nz, 1.0, 1.0, 1.0, (0.0,0.0,0.0), :uniform, "")
        
        grid = generate_grid(conf, dim)
        
        # Check sizes (including ghosts 2+2=4)
        @test grid.mx == nx + 4
        @test grid.my == ny + 4
        @test grid.mz == nz + 4
        
        # Check dx
        @test grid.dx ≈ 1.0/10
        @test grid.dy ≈ 1.0/10
        @test length(grid.z_face) == 10 + 5 # Nz+5
    end

    @testset "2. Memory/Fields" begin
        nx, ny, nz = 5, 5, 5
        est = estimate_memory_size(nx, ny, nz)
        @test est.bytes > 0
        @test est.arrays >= 20
        
        buffers = CFDBuffers(nx+4, ny+4, nz+4)
        @test size(buffers.u) == (nx+4, ny+4, nz+4)
    end

    @testset "3. WENO3 Unit Test" begin
        # 1D uniform test
        # Linear function u(x) = x should be reconstructed exactly by 3rd order
        # u_i = i * dx
        # Left reconstruction at i+1/2 (between u_i and u_i+1) should be (i+0.5)*dx
        
        dx = 0.1
        u_m1 = 1.0
        u_0  = 1.1
        u_p1 = 1.2
        
        # Ideal weights for smooth solution
        val = CleanroomSolver.Convection.weno3_reconstruct_left(u_m1, u_0, u_p1, dx, dx, dx, 1e-6)
        expected = 1.15
        
        @test isapprox(val, expected, atol=1e-5)
        
        # Shock case (Discontinuity)
        # u = 0, 0, 1. WENO should essentially pick 0 side.
        u_m1_s = 0.0
        u_0_s  = 0.0
        u_p1_s = 1.0
        
        val_s = CleanroomSolver.Convection.weno3_reconstruct_left(u_m1_s, u_0_s, u_p1_s, dx, dx, dx, 1e-6)
        # Check it is bounded and closer to 0 than linear interp 0.5
        @test val_s >= 0.0 && val_s < 0.5
    end
    
    @testset "4. Integration Test (Run Simulation)" begin
        # Path to test params
        param_file = joinpath(@__DIR__, "params.json")
        
        # Run simulation (10 steps)
        # Capture output or verify files exists
        
        # Clean up output
        out_dir = joinpath(@__DIR__, "output")
        rm(out_dir, recursive=true, force=true)
        
        run_simulation(param_file)
        
        @test isdir(out_dir)
        @test isfile(joinpath(out_dir, "history.txt"))
        
        # Verify monitoring data has 10 lines + header
        lines = readlines(joinpath(out_dir, "history.txt"))
        @test length(lines) >= 11 # Header + 10 steps
    end
end
