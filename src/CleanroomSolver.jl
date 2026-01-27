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
include("IO/Checkpoint.jl")
include("IO/SPHWriter.jl")
include("IO/InputReader.jl")
include("IO/Monitor.jl")

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

    nu_lam = 1.0 / dim_params.Re
    u0_nd = sim_params.initial_condition.velocity[1] / dim_params.U0
    v0_nd = sim_params.initial_condition.velocity[2] / dim_params.U0
    w0_nd = sim_params.initial_condition.velocity[3] / dim_params.U0
    U_ref = sqrt(u0_nd^2 + v0_nd^2 + w0_nd^2)
    dt_fixed = compute_dt(grid, sim_params.courant_number, nu_lam, U_ref)

    dx_min = min(grid.dx, grid.dy)
    if !isempty(grid.dz)
        dx_min = min(dx_min, minimum(grid.dz))
    end
    diff_num = nu_lam * dt_fixed / (dx_min^2)
    if diff_num >= 0.5
        println("Warning: Diffusion number D=$(diff_num) >= 0.5")
        if !sim_params.dry_run
            error("Diffusion number stability condition violated.")
        end
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
    bc_set = load_boundary_conditions(bc_path, dim_params)
    
    geo_path = replace(param_file, r"[^/\\]*$" => "geometry.json")
    if isfile(geo_path)
        objects = load_geometry(geo_path, dim_params)
        println("Loaded $(length(objects)) geometry objects.")
        fill_mask!(buffers.mask, objects, grid, "thread")
        apply_boundary_mask!(buffers.mask, grid, bc_set)
        apply_internal_boundary_mask!(buffers.mask, grid, bc_set)
        update_boundary_mask!(buffers.mask, grid, bc_set, 1.0)
    else
        println("No geometry file found. Domain is empty.")
        apply_boundary_mask!(buffers.mask, grid, bc_set)
        apply_internal_boundary_mask!(buffers.mask, grid, bc_set)
        update_boundary_mask!(buffers.mask, grid, bc_set, 1.0)
    end
    
    # Initialize fields (apply mask at solid regions)
    ic = sim_params.initial_condition
    u0, v0, w0 = ic.velocity ./ dim_params.U0
    p0 = ic.pressure / dim_params.U0^2
    
    initialize_fields!(buffers, grid, u0, v0, w0, p0)
    
    if isfile(geo_path)
        apply_object_velocity!(buffers.u, buffers.v, buffers.w, objects, grid)
    end
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
        sim_params.div_max_threshold,
        ndigits(sim_params.max_step) 
    )
    
    # Output directory handling
    # Create 'output' directory relative to the parameter file
    param_dir = dirname(abspath(param_file))
    out_dir = joinpath(param_dir, "output")
    mkpath(out_dir)
    println("Output Directory: $out_dir")
    
    init_monitor(
        monitor_config,
        joinpath(out_dir, "history.txt"),
        joinpath(out_dir, "condition.txt"),
        sim_params,
        grid,
        bc_set,
        dt_fixed,
        param_file
    )
    
    # Calculate fixed time step (dimensionless)
    println("Fixed time step (dimensionless): Δt* = $dt_fixed")
    
    # Initialize Ghost Cells with BCs before first step
    apply_velocity_cc_bcs!(
        buffers.u, buffers.v, buffers.w,
        grid, buffers.mask, bc_set, dt_fixed;
        reverse_flow_stabilization=sim_params.reverse_flow_stabilization
    )
    # Initialize face velocities from JSON initial velocity (u^n_face)
    apply_velocity_cf_bcs!(buffers.u_face_x, buffers.v_face_y, buffers.w_face_z, grid, buffers.mask, bc_set, dt_fixed)
    
    apply_pressure_bcs!(buffers.p, grid, buffers.mask, bc_set)
    BoundaryConditions.apply_internal_pressure_bcs!(buffers.p, grid, bc_set)

    println("Starting Time Loop...")
    step = start_step
    time = start_time
    
    # Allocate buffers for previous velocity (for dU calculation)
    u_prev = copy(buffers.u)
    v_prev = copy(buffers.v)
    w_prev = copy(buffers.w)
    
    t_conf = TimeConfig(time_scheme_enum, sim_params.courant_number, dt_fixed)

    cfl, cfl_i, cfl_j, cfl_k = 0.0, 0, 0, 0
    diff_num = 0.0
    u_max, u_i, u_j, u_k = 0.0, 0, 0, 0
    div_max, div_i, div_j, div_k = 0.0, 0, 0, 0
    
    while step < sim_params.max_step
        step += 1
        
        # Advance one time step
        pitr, pres = advance!(
            buffers, grid, bc_set, t_conf, sim_params.poisson, sim_params.smagorinsky_constant, nu_lam,
            sim_params.reverse_flow_stabilization, dt_fixed, "thread",
            rk_buffers=rk_buffers
        )
        time += dt_fixed

        # Stability monitoring (CFL/Diffusion) every step
        metrics = calculate_stability_metrics(buffers, grid, dt_fixed, nu_lam, sim_params.debug)
        cfl = metrics.cfl; diff_num = metrics.diff_num
        u_max = metrics.u_max; div_max = metrics.div_max
        div_i, div_j, div_k = metrics.div_max_loc
        
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
        
        mon_data = MonitorData(step, time, u_max, div_max, (div_i, div_j, div_k), dU, pitr, pres)
        open(joinpath(out_dir, "history.txt"), "a") do io
            log_step!(mon_data, monitor_config; console_io=stdout, history_io=io)
        end
        if sim_params.debug && monitor_config.console_interval > 0 && (step % monitor_config.console_interval == 0 || step == 1)
            log_flow_rates(buffers, grid, bc_set)
        end

        if sim_params.debug && monitor_config.console_interval > 0 && (step % monitor_config.console_interval == 0 || step == 1)
            uc = BoundaryConditions.compute_outflow_uc(
                buffers.u, buffers.v, buffers.w, grid, buffers.mask, bc_set
            )
            @printf("  Outflow Uc (nd): x_min=%.4e x_max=%.4e y_min=%.4e y_max=%.4e z_min=%.4e z_max=%.4e\n",
                uc.x_min, uc.x_max, uc.y_min, uc.y_max, uc.z_min, uc.z_max)
            if !isempty(bc_set.openings)
                for (op_idx, op) in enumerate(bc_set.openings)
                    if op.flow_type == OpeningOutlet && any(isnan, op.velocity)
                        face = op.boundary
                        region_check = (i, j, k) -> BoundaryConditions.is_in_opening(op, grid, i, j, k)
                        uc_region = BoundaryConditions.compute_region_average_velocity(
                            buffers.u, buffers.v, buffers.w, grid, face, region_check
                        )
                        label = isempty(op.name) ? "opening[$op_idx]" : op.name
                        @printf("  OpeningOutlet Uc (nd): %s @ %s = %.4e\n", label, face, uc_region)
                    end
                end
            end
        end

        if cfl > 1.0 || diff_num >= 0.5
            log_stability_violation(
                buffers, grid, bc_set, dt_fixed, metrics, sim_params.poisson, sim_params.debug
            )
            if !sim_params.dry_run
                break
            end
        end
        
        if check_divergence(div_max, monitor_config.div_threshold)
            # Convert to internal cell coordinates (excluding ghosts)
            internal_i = div_i - 2
            internal_j = div_j - 2
            internal_k = div_k - 2
            println("Error: Divergence detected at step $step (Div=$div_max)")
            println("  Max divergence location: internal cell ($internal_i, $internal_j, $internal_k)")
            println("  Grid indices (with ghost): ($div_i, $div_j, $div_k)")
            println("  Velocity at location: u=$(buffers.u[div_i,div_j,div_k]), v=$(buffers.v[div_i,div_j,div_k]), w=$(buffers.w[div_i,div_j,div_k])")
            break
        end

        # Time averaging (Welford) after start time
        if time >= sim_params.intervals.average_start_time_nd
            update_time_average!(buffers, "thread")
        end
        
        if sim_params.intervals.instantaneous > 0 && step % sim_params.intervals.instantaneous == 0
            fname = generate_sph_filename(joinpath(out_dir, "vel"), step)
            write_sph_vector(fname, buffers.u, buffers.v, buffers.w, grid, step, time, dim_params)
            fname_p = generate_sph_filename(joinpath(out_dir, "prs"), step)
            write_sph_scalar(fname_p, buffers.p, grid, step, time, dim_params)
        end
        
        if sim_params.intervals.checkpoint > 0 && step % sim_params.intervals.checkpoint == 0
            fname = generate_checkpoint_filename(step)
            write_checkpoint(joinpath(out_dir, fname), buffers, grid, step, time, dim_params)
        end
        
    end
    
    println("Simulation Completed.")
end



end # module CleanroomSolver
