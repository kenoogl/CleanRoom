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
    
    fill!(buffers.u, sim_params.initial_condition.velocity[1] / dim_params.U0) 
    fill!(buffers.v, sim_params.initial_condition.velocity[2] / dim_params.U0)
    fill!(buffers.w, sim_params.initial_condition.velocity[3] / dim_params.U0)
    fill!(buffers.p, sim_params.initial_condition.pressure / (dim_params.U0^2)) 
    
    # マスクを用いて固体領域の初期速度・圧力をゼロにする
    @inbounds @. buffers.u *= buffers.mask
    @inbounds @. buffers.v *= buffers.mask
    @inbounds @. buffers.w *= buffers.mask
    @inbounds @. buffers.p *= buffers.mask    
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
    apply_boundary_conditions!(buffers, grid, bc_set, dt_fixed, "thread")
    
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
            buffers, grid, bc_set, t_conf, sim_params.poisson, sim_params.smagorinsky_constant, nu_lam, dt_fixed, "thread",
            rk_buffers=rk_buffers
        )
        time += dt_fixed

        # Stability monitoring (CFL/Diffusion) every step
        cfl, cfl_i, cfl_j, cfl_k = compute_cfl(buffers, grid, dt_fixed)
        dx_min = min(grid.dx, grid.dy)
        if !isempty(grid.dz)
            dx_min = min(dx_min, minimum(grid.dz))
        end
        diff_num = nu_lam * dt_fixed / (dx_min^2)
        u_max, u_i, u_j, u_k = compute_u_max(buffers, grid)
        div_max, div_i, div_j, div_k = compute_divergence_max(buffers, grid)
        
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
            # Inflow/Outflow face flow rates (outward normal positive)
            flow_xmin = 0.0
            flow_xmax = 0.0
            flow_ymin = 0.0
            flow_ymax = 0.0
            flow_zmin = 0.0
            flow_zmax = 0.0
            flow_xmin_in = 0.0
            flow_xmax_in = 0.0
            flow_ymin_in = 0.0
            flow_ymax_in = 0.0
            flow_zmin_in = 0.0
            flow_zmax_in = 0.0
            flow_xmin_out = 0.0
            flow_xmax_out = 0.0
            flow_ymin_out = 0.0
            flow_ymax_out = 0.0
            flow_zmin_out = 0.0
            flow_zmax_out = 0.0
            un_xmin_min = Inf
            un_xmin_max = -Inf
            un_xmax_min = Inf
            un_xmax_max = -Inf
            un_ymin_min = Inf
            un_ymin_max = -Inf
            un_ymax_min = Inf
            un_ymax_max = -Inf
            un_zmin_min = Inf
            un_zmin_max = -Inf
            un_zmax_min = Inf
            un_zmax_max = -Inf
            x_min_io = bc_set.x_min.velocity_type == Inflow || bc_set.x_min.velocity_type == Outflow
            x_max_io = bc_set.x_max.velocity_type == Inflow || bc_set.x_max.velocity_type == Outflow
            y_min_io = bc_set.y_min.velocity_type == Inflow || bc_set.y_min.velocity_type == Outflow
            y_max_io = bc_set.y_max.velocity_type == Inflow || bc_set.y_max.velocity_type == Outflow
            z_min_io = bc_set.z_min.velocity_type == Inflow || bc_set.z_min.velocity_type == Outflow
            z_max_io = bc_set.z_max.velocity_type == Inflow || bc_set.z_max.velocity_type == Outflow
            if x_min_io
                @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2
                    un = -buffers.u_face_x[3, j, k]  # outward normal for x_min
                    area = grid.dy * grid.dz[k]
                    flow_xmin += un * area
                    if un >= 0.0
                        flow_xmin_out += un * area
                    else
                        flow_xmin_in += -un * area
                    end
                    un_xmin_min = min(un_xmin_min, un)
                    un_xmin_max = max(un_xmin_max, un)
                end
            end
            if x_max_io
                @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2
                    un = buffers.u_face_x[grid.mx-1, j, k]  # outward normal for x_max
                    area = grid.dy * grid.dz[k]
                    flow_xmax += un * area
                    if un >= 0.0
                        flow_xmax_out += un * area
                    else
                        flow_xmax_in += -un * area
                    end
                    un_xmax_min = min(un_xmax_min, un)
                    un_xmax_max = max(un_xmax_max, un)
                end
            end
            if y_min_io
                @inbounds for k in 3:grid.mz-2, i in 3:grid.mx-2
                    un = -buffers.v_face_y[i, 3, k]  # outward normal for y_min
                    area = grid.dx * grid.dz[k]
                    flow_ymin += un * area
                    if un >= 0.0
                        flow_ymin_out += un * area
                    else
                        flow_ymin_in += -un * area
                    end
                    un_ymin_min = min(un_ymin_min, un)
                    un_ymin_max = max(un_ymin_max, un)
                end
            end
            if y_max_io
                @inbounds for k in 3:grid.mz-2, i in 3:grid.mx-2
                    un = buffers.v_face_y[i, grid.my-1, k]  # outward normal for y_max
                    area = grid.dx * grid.dz[k]
                    flow_ymax += un * area
                    if un >= 0.0
                        flow_ymax_out += un * area
                    else
                        flow_ymax_in += -un * area
                    end
                    un_ymax_min = min(un_ymax_min, un)
                    un_ymax_max = max(un_ymax_max, un)
                end
            end
            if z_min_io
                @inbounds for j in 3:grid.my-2, i in 3:grid.mx-2
                    un = -buffers.w_face_z[i, j, 3]  # outward normal for z_min
                    area = grid.dx * grid.dy
                    flow_zmin += un * area
                    if un >= 0.0
                        flow_zmin_out += un * area
                    else
                        flow_zmin_in += -un * area
                    end
                    un_zmin_min = min(un_zmin_min, un)
                    un_zmin_max = max(un_zmin_max, un)
                end
            end
            if z_max_io
                @inbounds for j in 3:grid.my-2, i in 3:grid.mx-2
                    un = buffers.w_face_z[i, j, grid.mz-1]  # outward normal for z_max
                    area = grid.dx * grid.dy
                    flow_zmax += un * area
                    if un >= 0.0
                        flow_zmax_out += un * area
                    else
                        flow_zmax_in += -un * area
                    end
                    un_zmax_min = min(un_zmax_min, un)
                    un_zmax_max = max(un_zmax_max, un)
                end
            end
            xmin_str = x_min_io ? @sprintf("%.4e", flow_xmin) : "N/A"
            xmax_str = x_max_io ? @sprintf("%.4e", flow_xmax) : "N/A"
            ymin_str = y_min_io ? @sprintf("%.4e", flow_ymin) : "N/A"
            ymax_str = y_max_io ? @sprintf("%.4e", flow_ymax) : "N/A"
            zmin_str = z_min_io ? @sprintf("%.4e", flow_zmin) : "N/A"
            zmax_str = z_max_io ? @sprintf("%.4e", flow_zmax) : "N/A"
            println("  Flow (nd, outward+): x_min=$xmin_str x_max=$xmax_str y_min=$ymin_str y_max=$ymax_str z_min=$zmin_str z_max=$zmax_str")
            println("  FaceFlow split (nd):",
                    x_min_io ? @sprintf(" x_min[out=%.4e in=%.4e un(%.3e..%.3e)]", flow_xmin_out, flow_xmin_in, un_xmin_min, un_xmin_max) : " x_min[N/A]",
                    x_max_io ? @sprintf(" x_max[out=%.4e in=%.4e un(%.3e..%.3e)]", flow_xmax_out, flow_xmax_in, un_xmax_min, un_xmax_max) : " x_max[N/A]",
                    y_min_io ? @sprintf(" y_min[out=%.4e in=%.4e un(%.3e..%.3e)]", flow_ymin_out, flow_ymin_in, un_ymin_min, un_ymin_max) : " y_min[N/A]",
                    y_max_io ? @sprintf(" y_max[out=%.4e in=%.4e un(%.3e..%.3e)]", flow_ymax_out, flow_ymax_in, un_ymax_min, un_ymax_max) : " y_max[N/A]",
                    z_min_io ? @sprintf(" z_min[out=%.4e in=%.4e un(%.3e..%.3e)]", flow_zmin_out, flow_zmin_in, un_zmin_min, un_zmin_max) : " z_min[N/A]",
                    z_max_io ? @sprintf(" z_max[out=%.4e in=%.4e un(%.3e..%.3e)]", flow_zmax_out, flow_zmax_in, un_zmax_min, un_zmax_max) : " z_max[N/A]")
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
            if sim_params.debug
                u_c = buffers.u[cfl_i, cfl_j, cfl_k]
                v_c = buffers.v[cfl_i, cfl_j, cfl_k]
                w_c = buffers.w[cfl_i, cfl_j, cfl_k]
                m_c = buffers.mask[cfl_i, cfl_j, cfl_k]
                u_m = buffers.u[u_i, u_j, u_k]
                v_m = buffers.v[u_i, u_j, u_k]
                w_m = buffers.w[u_i, u_j, u_k]
                ustar_max = 0.0
                us_i, us_j, us_k = 0, 0, 0
                rhs_max = 0.0
                rh_i, rh_j, rh_k = 0, 0, 0
                p_max = 0.0
                p_i, p_j, p_k = 0, 0, 0
                divu_max = 0.0
                du_i, du_j, du_k = 0, 0, 0
                @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2, i in 3:grid.mx-2
                    m0 = buffers.mask[i, j, k]
                    if m0 == 0.0
                        continue
                    end
                    um = buffers.u_star[i, j, k]
                    vm = buffers.v_star[i, j, k]
                    wm = buffers.w_star[i, j, k]
                    mag = sqrt(um^2 + vm^2 + wm^2)
                    if mag > ustar_max
                        ustar_max = mag
                        us_i, us_j, us_k = i, j, k
                    end
                    rv = abs(buffers.rhs[i, j, k])
                    if rv > rhs_max
                        rhs_max = rv
                        rh_i, rh_j, rh_k = i, j, k
                    end
                    pv = abs(buffers.p[i, j, k])
                    if pv > p_max
                        p_max = pv
                        p_i, p_j, p_k = i, j, k
                    end
                    # div(u*) max (use face-averaged u* with mask)
                    u_fr = 0.5 * (buffers.u_star[i+1, j, k] + buffers.u_star[i, j, k]) * buffers.mask[i+1, j, k] * m0
                    u_fl = 0.5 * (buffers.u_star[i, j, k] + buffers.u_star[i-1, j, k]) * buffers.mask[i-1, j, k] * m0
                    v_ft = 0.5 * (buffers.v_star[i, j+1, k] + buffers.v_star[i, j, k]) * buffers.mask[i, j+1, k] * m0
                    v_fb = 0.5 * (buffers.v_star[i, j, k] + buffers.v_star[i, j-1, k]) * buffers.mask[i, j-1, k] * m0
                    w_fu = 0.5 * (buffers.w_star[i, j, k+1] + buffers.w_star[i, j, k]) * buffers.mask[i, j, k+1] * m0
                    w_fd = 0.5 * (buffers.w_star[i, j, k] + buffers.w_star[i, j, k-1]) * buffers.mask[i, j, k-1] * m0
                    dz = grid.dz[k]
                    divu = (u_fr - u_fl) / grid.dx + (v_ft - v_fb) / grid.dy + (w_fu - w_fd) / dz
                    abs_divu = abs(divu)
                    if abs_divu > divu_max
                        divu_max = abs_divu
                        du_i, du_j, du_k = i, j, k
                    end
                end
                # Poisson residuals using pressure mask (inflow/outflow/opening set to 0)
                update_boundary_mask!(buffers.mask, grid, bc_set, 0.0)
                alpha = sim_params.poisson.mach2 / (dt_fixed * dt_fixed)
                res0 = PressureSolver.compute_residual_sor(buffers.p_prev, buffers.rhs, buffers.mask, grid, sim_params.poisson.omega, alpha)
                res1 = PressureSolver.compute_residual_sor(buffers.p, buffers.rhs, buffers.mask, grid, sim_params.poisson.omega, alpha)
                update_boundary_mask!(buffers.mask, grid, bc_set, 1.0)
                println("Stability violation: CFL=$cfl, D=$diff_num")
                println("  CFL max at (i,j,k)=($(cfl_i),$(cfl_j),$(cfl_k)) mask=$(m_c) u=$(u_c) v=$(v_c) w=$(w_c)")
                println("  Umax at (i,j,k)=($(u_i),$(u_j),$(u_k)) |U|=$(u_max) u=$(u_m) v=$(v_m) w=$(w_m)")
                println("  U* max at (i,j,k)=($(us_i),$(us_j),$(us_k)) |U*|=$(ustar_max)")
                println("  RHS max at (i,j,k)=($(rh_i),$(rh_j),$(rh_k)) |rhs|=$(rhs_max)")
                println("  |p| max at (i,j,k)=($(p_i),$(p_j),$(p_k)) |p|=$(p_max)")
                println("  |div(u*)| max at (i,j,k)=($(du_i),$(du_j),$(du_k)) |div|=$(divu_max)")
                println("  Poisson residuals: r0=$(res0) r1=$(res1)")
                println("  Poisson: iter=$(pitr) res=$(pres)")
            else
                println("Stability violation: CFL=$cfl, D=$diff_num")
            end
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
