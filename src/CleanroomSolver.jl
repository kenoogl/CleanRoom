module CleanroomSolver

# --- Includes (Order matters) ---
include("Common.jl")
include("Grid.jl")
include("Fields.jl")
include("BC/BoundaryConditions.jl")
include("BC/Geometry.jl")
include("Physics/PressureSolver.jl")
include("Physics/Turbulence.jl")
include("Physics/Diffusion.jl")
include("Physics/Convection.jl")
include("Physics/FractionalStep.jl")
include("TimeIntegration.jl") 
include("IO/Visualization.jl")
include("IO/Monitor.jl")
include("IO/Checkpoint.jl")
include("IO/SPHWriter.jl")
include("IO/InputReader.jl")

# --- Using submodules ---
using .Common
using .Grid
using .Fields
using .BoundaryConditions
using .Geometry
using .PressureSolver
using .Turbulence
using .Diffusion
using .Convection
using .FractionalStep
using .TimeIntegration
using .Visualization
using .Monitor
using .Checkpoint
using .SPHWriter
using .InputReader

using Printf

export run_simulation

"""
    run_simulation(param_file::String)

シミュレーションを実行するメインエントリーポイント。
"""
function run_simulation(param_file::String)
    println("--- CleanroomSolver Started ---")
    println("Loading parameters from: $param_file")
    
    sim_params = load_parameters(param_file)
    dim_params = sim_params.dim_params
    
    println("Parameters loaded. L0=$(dim_params.L0), U0=$(dim_params.U0), Re=$(dim_params.Re)")
    
    println("Generating Grid...")
    grid = generate_grid(sim_params.grid_config, dim_params)
    println("Grid: $(grid.mx)x$(grid.my)x$(grid.mz) (including ghosts)")
    
    mem_ok, mem_msg = check_memory_availability(sim_params.grid_config.Nx, sim_params.grid_config.Ny, sim_params.grid_config.Nz)
    println("Memory Check: $mem_msg")
    if !mem_ok && !sim_params.dry_run
        error("Insufficient memory.")
    end
    
    if sim_params.dry_run
        println("Dry run completed.")
        return
    end
    
    println("Allocating Fields...")
    buffers = CFDBuffers(grid.mx, grid.my, grid.mz)
    
    # RK buffering
    time_scheme_enum = if sim_params.time_scheme == :RK2
        RK2
    elseif sim_params.time_scheme == :RK4
        RK4
    else
        Euler
    end
    
    rk_buffers = allocate_rk_buffers(Symbol(sim_params.time_scheme), grid.mx, grid.my, grid.mz)
    if !isnothing(rk_buffers)
        println("Allocated RK Buffers for $(sim_params.time_scheme)")
    end
    
    bc_path = replace(param_file, r"[^/\\]*$" => "boundary_conditions.json")
    if !isfile(bc_path)
         println("Warning: BC file not found at $bc_path.")
         return
    end
    bc_set = load_boundary_conditions(bc_path)
    
    geo_path = replace(param_file, r"[^/\\]*$" => "geometry.json")
    if isfile(geo_path)
        objects = load_geometry(geo_path)
        println("Loaded $(length(objects)) geometry objects.")
        fill_mask!(buffers.mask, objects, grid, "thread")
    else
        println("No geometry file found. Domain is empty.")
    end
    
    fill!(buffers.u, sim_params.initial_condition.velocity[1] / dim_params.U0) 
    fill!(buffers.v, sim_params.initial_condition.velocity[2] / dim_params.U0)
    fill!(buffers.w, sim_params.initial_condition.velocity[3] / dim_params.U0)
    fill!(buffers.p, sim_params.initial_condition.pressure / (dim_params.U0^2)) 
    
    start_step = 0
    start_time = 0.0
    if sim_params.start == :restart
        println("Restarting from $(sim_params.restart_file)...")
        start_step, start_time = read_checkpoint(sim_params.restart_file, buffers, dim_params)
        println("Restarted at Step $start_step, Time $start_time")
    end
    
    monitor_config = MonitorConfig(
        sim_params.intervals.display,
        sim_params.intervals.history,
        1e5,
        ndigits(sim_params.max_step) 
    )
    
    # Output directory handling
    # Create 'output' directory relative to the parameter file
    param_dir = dirname(abspath(param_file))
    out_dir = joinpath(param_dir, "output")
    mkpath(out_dir)
    println("Output Directory: $out_dir")
    
    init_monitor(monitor_config, joinpath(out_dir, "history.txt"))
    
    # Calculate fixed time step (dimensionless)
    nu_lam = 1.0 / dim_params.Re
    dt_fixed = compute_dt(buffers, grid, sim_params.courant_number, nu_lam)
    println("Fixed time step (dimensionless): Δt* = $dt_fixed")
    
    # Output calculation conditions to condition.txt
    write_condition_file(out_dir, param_file, dim_params, sim_params, grid, dt_fixed)
    
    println("Starting Time Loop...")
    step = start_step
    time = start_time
    
    # Allocate buffers for previous velocity (for dU calculation)
    u_prev = copy(buffers.u)
    v_prev = copy(buffers.v)
    w_prev = copy(buffers.w)
    
    t_conf = TimeConfig(time_scheme_enum, sim_params.courant_number, dt_fixed)
    
    while step < sim_params.max_step
        step += 1
        
        # Advance one time step
        pitr, pres = advance!(
            buffers, grid, bc_set, t_conf, sim_params.poisson, 0.1, nu_lam, dt_fixed, "thread",
            rk_buffers=rk_buffers
        )
        time += dt_fixed
        
        u_max = compute_u_max(buffers, grid)
        div = compute_divergence_max(buffers, grid)
        
        # Compute velocity change for steady-state check
        # dU = sqrt( Σ( (u^n+1 - u^n)² + (v^n+1 - v^n)² + (w^n+1 - w^n)² ) )
        dU_sum = 0.0
        @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2, i in 3:grid.mx-2
            dU_sum += (buffers.u[i,j,k] - u_prev[i,j,k])^2 +
                      (buffers.v[i,j,k] - v_prev[i,j,k])^2 +
                      (buffers.w[i,j,k] - w_prev[i,j,k])^2
        end
        dU = sqrt(dU_sum)
        
        # Store current velocity as previous for next step
        u_prev .= buffers.u
        v_prev .= buffers.v
        w_prev .= buffers.w
        
        mon_data = MonitorData(step, time, u_max, div, dU, pitr, pres)
        log_step!(mon_data, monitor_config, stdout)
        open(joinpath(out_dir, "history.txt"), "a") do io
            log_step!(mon_data, monitor_config, io)
        end
        
        if check_divergence(div, monitor_config.div_threshold)
            println("Error: Divergence detected at step $step (Div=$div)")
            break
        end 
        
        if sim_params.intervals.instantaneous > 0 && step % sim_params.intervals.instantaneous == 0
            fname = generate_sph_filename(joinpath(out_dir, "vel"), step)
            write_sph_vector(fname, buffers.u, buffers.v, buffers.w, grid, step, time, sim_params.output_dimensional, dim_params)
            fname_p = generate_sph_filename(joinpath(out_dir, "prs"), step)
            write_sph_scalar(fname_p, buffers.p, grid, step, time, sim_params.output_dimensional, dim_params)
        end
        
        if sim_params.intervals.checkpoint > 0 && step % sim_params.intervals.checkpoint == 0
            fname = generate_checkpoint_filename(step)
            write_checkpoint(joinpath(out_dir, fname), buffers, grid, step, time, false, dim_params)
        end
        
        if sim_params.visualization.interval > 0 && step % sim_params.visualization.interval == 0
             render_slice(buffers, grid, sim_params.visualization, step, dim_params)
        end
        
    end
    
    println("Simulation Completed.")
end


"""
    write_condition_file(out_dir, param_file, dim_params, sim_params, grid, dt_fixed)

計算条件をファイルに出力する。
"""
function write_condition_file(out_dir, param_file, dim_params, sim_params, grid, dt_fixed)
    condition_file = joinpath(out_dir, "condition.txt")
    open(condition_file, "w") do io
        println(io, "=== Calculation Conditions ===")
        println(io, "Parameter File: $(abspath(param_file))")
        println(io, "")
        println(io, "--- Physical Parameters ---")
        @printf(io, "  %-22s %12.4g [m]    \n", "Reference Length L0:", dim_params.L0)
        @printf(io, "  %-22s %12.4g [m/s]  \n", "Reference Velocity U0:", dim_params.U0)
        @printf(io, "  %-22s %12.4e [m²/s] \n", "Kinematic Viscosity ν:", dim_params.nu)
        @printf(io, "  %-22s %12.4g        \n", "Reynolds Number Re:", dim_params.Re)
        @printf(io, "  %-22s %12.4g [s]    \n", "Reference Time T0:", dim_params.T0)
        println(io, "")
        println(io, "--- Grid Parameters ---")
        @printf(io, "  %-22s (%d, %d, %d)\n", "Grid Size (Nx,Ny,Nz):", grid.mx-4, grid.my-4, grid.mz-4)
        @printf(io, "  %-22s %12.4g [m]\n", "Domain Lx:", sim_params.grid_config.Lx)
        @printf(io, "  %-22s %12.4g [m]\n", "Domain Ly:", sim_params.grid_config.Ly)
        Lz = (grid.z_face[end-2] - grid.z_face[3]) * dim_params.L0
        @printf(io, "  %-22s %12.4g [m]\n", "Domain Lz:", Lz)
        @printf(io, "  %-22s %12.6g [m]\n", "Cell Δx:", grid.dx * dim_params.L0)
        @printf(io, "  %-22s %12.6g [m]\n", "Cell Δy:", grid.dy * dim_params.L0)
        dz_min = minimum(grid.dz[3:end-2]) * dim_params.L0
        dz_max = maximum(grid.dz[3:end-2]) * dim_params.L0
        @printf(io, "  %-22s %12.6g [m] (min)\n", "Cell Δz:", dz_min)
        @printf(io, "  %-22s %12.6g [m] (max)\n", "", dz_max)
        println(io, "")
        println(io, "--- Time Integration ---")
        @printf(io, "  %-22s %s\n", "Time Scheme:", sim_params.time_scheme)
        @printf(io, "  %-22s %12.4g\n", "Courant Number:", sim_params.courant_number)
        @printf(io, "  %-22s %12.6g (dimensionless)\n", "Fixed Time Step Δt*:", dt_fixed)
        @printf(io, "  %-22s %12.6g [s]\n", "Fixed Time Step Δt:", dt_fixed * dim_params.T0)
        @printf(io, "  %-22s %12d\n", "Max Steps:", sim_params.max_step)
        total_time_nd = dt_fixed * sim_params.max_step
        total_time_dim = total_time_nd * dim_params.T0
        @printf(io, "  %-22s %12.4g [s]\n", "Total Simulation Time:", total_time_dim)
        println(io, "")
        println(io, "--- Poisson Solver ---")
        @printf(io, "  %-22s %s\n", "Solver:", sim_params.poisson.solver)
        @printf(io, "  %-22s %12.1e\n", "Convergence Criteria:", sim_params.poisson.tol)
        @printf(io, "  %-22s %12d\n", "Max Iterations:", sim_params.poisson.max_iter)
        println(io, "")
        println(io, "--- Output Intervals ---")
        @printf(io, "  %-22s %8d steps\n", "Display:", sim_params.intervals.display)
        @printf(io, "  %-22s %8d steps\n", "History:", sim_params.intervals.history)
        @printf(io, "  %-22s %8d steps\n", "Instantaneous:", sim_params.intervals.instantaneous)
        @printf(io, "  %-22s %8d steps\n", "Checkpoint:", sim_params.intervals.checkpoint)
    end
    println("Condition file: $condition_file")
end

end # module CleanroomSolver
