module Monitor

using Printf
using ..Common
using ..Fields
using ..Grid
using ..InputReader
using ..BoundaryConditions
using ..PressureSolver

export MonitorData, MonitorConfig, init_monitor, log_step!, check_divergence
export compute_u_max, compute_cfl, compute_divergence_max
export calculate_stability_metrics, log_flow_rates, log_stability_violation

struct MonitorData
    step::Int
    time::Float64
    u_max::Float64
    div_max::Float64
    div_loc::NTuple{3, Int}  # (i, j, k) of max divergence
    dU::Float64
    pressure_itr::Int
    pressure_residual::Float64
end
struct MonitorConfig
    console_interval::Int
    history_interval::Int
    div_threshold::Float64
    step_digits::Int
end

function init_monitor(
    config::MonitorConfig,
    history_path::String,
    condition_path::String,
    sim_params::SimulationParams,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    dt_fixed::Float64,
    param_file::String
)
    # Write condition file
    write_condition_file(condition_path, param_file, sim_params, grid, bc_set, dt_fixed)

    open(history_path, "w") do io
        step_width = max(config.step_digits, 4)
        s_step = rpad("step", step_width)
        
        print(io, s_step, " ")
        # Header with Loc column
        @printf(io, "%14s %12s %12s %12s %12s %5s %13s\n",
            "time", "Umax", "divMax", "Loc(i,j,k)", "dU", "ItrP", "ResP")
    end
end

function write_condition_file(
    condition_path::String,
    param_file::String,
    sim_params::SimulationParams,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    dt_fixed::Float64
)
    dim_params = sim_params.dim_params
    open(condition_path, "w") do io
        println(io, "=== Calculation Conditions ===")
        println(io, "Parameter File: $(abspath(param_file))")
        println(io, "")
        println(io, "--- Physical Parameters ---")
        @printf(io, "  %-22s %12.4g [m]    \n", "Reference Length L0:", dim_params.L0)
        @printf(io, "  %-22s %12.4g [m/s]  \n", "Reference Velocity U0:", dim_params.U0)
        @printf(io, "  %-22s %12.4e [m²/s] \n", "Kinematic Viscosity ν:", dim_params.nu)
        @printf(io, "  %-22s %12.4g        \n", "Reynolds Number Re:", dim_params.Re)
        @printf(io, "  %-22s %12.4g [s]    \n", "Reference Time T0:", dim_params.T0)
        @printf(io, "  %-22s %12.4g        \n", "Mach Number M:", sqrt(sim_params.poisson.mach2))
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
        @printf(io, "  %-22s %s\n", "Debug:", sim_params.debug ? "yes" : "no")
        @printf(io, "  %-22s %s\n", "Reverse Flow Stabilization:", sim_params.reverse_flow_stabilization ? "yes" : "no")
        
        # Diffusion number: (1/Re) * dt / dx^2
        inv_re = 1.0 / dim_params.Re
        diff_x = inv_re * dt_fixed / (grid.dx^2)
        diff_y = inv_re * dt_fixed / (grid.dy^2)
        dz_min = minimum(grid.dz[3:end-2])
        diff_z = inv_re * dt_fixed / (dz_min^2)
        diff_max = max(diff_x, diff_y, diff_z)
        @printf(io, "  %-22s %12.4g\n", "Diffusion Number:", diff_max)
        
        @printf(io, "  %-22s %12d\n", "Max Steps:", sim_params.max_step)
        total_time_nd = dt_fixed * sim_params.max_step
        total_time_dim = total_time_nd * dim_params.T0
        @printf(io, "  %-22s %12.4g [s]\n", "Total Simulation Time:", total_time_dim)
        println(io, "")
        
        println(io, "--- Boundary Conditions ---")
        function print_ext(label, bc)
            val_str = ""
            if bc.velocity_type == Inflow || bc.velocity_type == SlidingWall
                v = bc.velocity_value
                vx, vy, vz = v .* dim_params.U0
                val_str = @sprintf("(%.2f, %.2f, %.2f)", vx, vy, vz)
            end
            @printf(io, "  %-10s %-12s %s\n", label, bc.velocity_type, val_str)
        end
        print_ext("x_min:", bc_set.x_min)
        print_ext("x_max:", bc_set.x_max)
        print_ext("y_min:", bc_set.y_min)
        print_ext("y_max:", bc_set.y_max)
        print_ext("z_min:", bc_set.z_min)
        print_ext("z_max:", bc_set.z_max)
        @printf(io, "  %-22s %d\n", "Openings:", length(bc_set.openings))
        @printf(io, "  %-22s %d\n", "Internal Boundaries:", length(bc_set.internal_boundaries))
        println(io, "")

        println(io, "--- Poisson Solver ---")
        @printf(io, "  %-22s %s\n", "Solver:", sim_params.poisson.solver)
        if sim_params.poisson.solver == CG || sim_params.poisson.solver == BiCGSTAB
            @printf(io, "  %-22s %s\n", "Preconditioner:", sim_params.poisson.preconditioner)
        else
            @printf(io, "  %-22s %s\n", "Preconditioner:", "N/A")
        end
        if sim_params.poisson.solver == RedBlackSOR
            @printf(io, "  %-22s %12.4g\n", "SOR Omega:", sim_params.poisson.omega)
        end
        @printf(io, "  %-22s %12.1e\n", "Convergence Criteria:", sim_params.poisson.tol)
        @printf(io, "  %-22s %12d\n", "Max Iterations:", sim_params.poisson.max_iter)
        println(io, "")
        println(io, "--- Output Intervals ---")
        @printf(io, "  %-22s %8d steps\n", "Display:", sim_params.intervals.display)
        @printf(io, "  %-22s %8d steps\n", "History:", sim_params.intervals.history)
        @printf(io, "  %-22s %8d steps\n", "Instantaneous:", sim_params.intervals.instantaneous)
        @printf(io, "  %-22s %8d steps\n", "Checkpoint:", sim_params.intervals.checkpoint)
    end
    println("Condition file: $condition_path")
end

function compute_u_max(buffers::CFDBuffers, grid::GridData)
    umax = 0.0
    u_i, u_j, u_k = 0, 0, 0
    u, v, w = buffers.u, buffers.v, buffers.w
    @inbounds for k in 1:grid.mz, j in 1:grid.my, i in 1:grid.mx
        mag = sqrt(u[i, j, k]^2 + v[i, j, k]^2 + w[i, j, k]^2)
        if mag > umax
            umax = mag
            u_i, u_j, u_k = i, j, k
        end
    end
    return umax, u_i, u_j, u_k
end

function compute_cfl(buffers::CFDBuffers, grid::GridData, dt::Float64)
    cfl_max = 0.0
    cfl_i, cfl_j, cfl_k = 0, 0, 0
    u, v, w = buffers.u, buffers.v, buffers.w
    mask = buffers.mask
    dx, dy = grid.dx, grid.dy
    @inbounds for k in 3:grid.mz-2
        dz = grid.dz[k]
        for j in 3:grid.my-2, i in 3:grid.mx-2
             m0 = mask[i, j, k]
             cfl_x = abs(u[i, j, k]) * dt / dx
             cfl_y = abs(v[i, j, k]) * dt / dy
             cfl_z = abs(w[i, j, k]) * dt / dz
             val = max(cfl_x, cfl_y, cfl_z) * m0
             if val > cfl_max
                 cfl_max = val
                 cfl_i, cfl_j, cfl_k = i, j, k
             end
        end
    end
    return cfl_max, cfl_i, cfl_j, cfl_k
end

function compute_divergence_max(buffers::CFDBuffers, grid::GridData)::Tuple{Float64, Int, Int, Int}
    max_div = 0.0
    max_i, max_j, max_k = 0, 0, 0
    
    u_face = buffers.u_face_x
    v_face = buffers.v_face_y
    w_face = buffers.w_face_z
    
    @inbounds for k in 3:grid.mz-2
        dz = grid.z_face[k+1] - grid.z_face[k]
        for j in 3:grid.my-2, i in 3:grid.mx-2
            div = (u_face[i+1, j, k] - u_face[i, j, k]) / grid.dx +
                  (v_face[i, j+1, k] - v_face[i, j, k]) / grid.dy +
                  (w_face[i, j, k+1] - w_face[i, j, k]) / dz
            
            abs_div = abs(div)
            if abs_div > max_div
                max_div = abs_div
                max_i, max_j, max_k = i, j, k
            end
        end
    end
    return (max_div, max_i, max_j, max_k)
end

function log_step!(
    data::MonitorData,
    config::MonitorConfig;
    console_io::Union{IO, Nothing}=stdout,
    history_io::Union{IO, Nothing}=nothing
)
    step_width = max(config.step_digits, 4)
    
    # Format location string: "(i,j,k)"
    s_loc = @sprintf("(%d,%d,%d)", data.div_loc[1], data.div_loc[2], data.div_loc[3])

    # Console output
    if console_io !== nothing && (data.step % config.console_interval == 0 || data.step == 1)
        if data.step == 1
            @printf(console_io, "%*s %14s %12s %12s %12s %12s %5s %13s\n",
                step_width, "Step", "Time", "Umax", "Div", "Loc", "dU", "ItrP", "ResP")
        end
        @printf(console_io, "%*d %14.6e %12.4e %12.4e %12s %12.4e %5d %13.5e\n",
            step_width, data.step, data.time, data.u_max, data.div_max, s_loc,
            data.dU, data.pressure_itr, data.pressure_residual)
    end

    if history_io !== nothing && (data.step % config.history_interval == 0 || data.step == 1)
        s_step = rpad(string(data.step), step_width)
        print(history_io, s_step, " ")
        
        @printf(history_io, "%14.6e %12.4e %12.4e %12s %12.4e %5d %13.5e\n",
            data.time, data.u_max, data.div_max, s_loc, data.dU, data.pressure_itr, data.pressure_residual)
        flush(history_io)
    end
end

function check_divergence(div_max::Float64, threshold::Float64)::Bool
    return div_max > threshold || isnan(div_max) || isinf(div_max)
end

function calculate_stability_metrics(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    nu_lam::Float64,
    debug::Bool
)
    if debug
        cfl, cfl_i, cfl_j, cfl_k = compute_cfl(buffers, grid, dt)
        dx_min = min(grid.dx, grid.dy)
        if !isempty(grid.dz)
             dx_min = min(dx_min, minimum(grid.dz))
        end
        diff_num = nu_lam * dt / (dx_min^2)
        u_max, u_i, u_j, u_k = compute_u_max(buffers, grid)
        div_max, div_i, div_j, div_k = compute_divergence_max(buffers, grid)
        
        return (
            cfl=cfl, cfl_loc=(cfl_i, cfl_j, cfl_k),
            diff_num=diff_num,
            u_max=u_max, u_max_loc=(u_i, u_j, u_k),
            div_max=div_max, div_max_loc=(div_i, div_j, div_k)
        )
    else
        return (
            cfl=0.0, cfl_loc=(0,0,0),
            diff_num=0.0,
            u_max=0.0, u_max_loc=(0,0,0),
            div_max=0.0, div_max_loc=(0,0,0)
        )
    end
end

function log_flow_rates(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet
)
    # Inflow/Outflow face flow rates (outward normal positive)
    flow_xmin = 0.0; flow_xmax = 0.0
    flow_ymin = 0.0; flow_ymax = 0.0
    flow_zmin = 0.0; flow_zmax = 0.0
    
    flow_xmin_in = 0.0; flow_xmax_in = 0.0
    flow_ymin_in = 0.0; flow_ymax_in = 0.0
    flow_zmin_in = 0.0; flow_zmax_in = 0.0
    
    flow_xmin_out = 0.0; flow_xmax_out = 0.0
    flow_ymin_out = 0.0; flow_ymax_out = 0.0
    flow_zmin_out = 0.0; flow_zmax_out = 0.0
    
    un_xmin_min = Inf; un_xmin_max = -Inf
    un_xmax_min = Inf; un_xmax_max = -Inf
    un_ymin_min = Inf; un_ymin_max = -Inf
    un_ymax_min = Inf; un_ymax_max = -Inf
    un_zmin_min = Inf; un_zmin_max = -Inf
    un_zmax_min = Inf; un_zmax_max = -Inf
    
    x_min_io = bc_set.x_min.velocity_type == Inflow || bc_set.x_min.velocity_type == Outflow
    x_max_io = bc_set.x_max.velocity_type == Inflow || bc_set.x_max.velocity_type == Outflow
    y_min_io = bc_set.y_min.velocity_type == Inflow || bc_set.y_min.velocity_type == Outflow
    y_max_io = bc_set.y_max.velocity_type == Inflow || bc_set.y_max.velocity_type == Outflow
    z_min_io = bc_set.z_min.velocity_type == Inflow || bc_set.z_min.velocity_type == Outflow
    z_max_io = bc_set.z_max.velocity_type == Inflow || bc_set.z_max.velocity_type == Outflow
    
    if x_min_io
        @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2
            un = -buffers.u_face_x[3, j, k]
            area = grid.dy * grid.dz[k]
            flow_xmin += un * area
            if un >= 0.0; flow_xmin_out += un * area; else; flow_xmin_in += -un * area; end
            un_xmin_min = min(un_xmin_min, un); un_xmin_max = max(un_xmin_max, un)
        end
    end
    if x_max_io
        @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2
            un = buffers.u_face_x[grid.mx-1, j, k]
            area = grid.dy * grid.dz[k]
            flow_xmax += un * area
            if un >= 0.0; flow_xmax_out += un * area; else; flow_xmax_in += -un * area; end
            un_xmax_min = min(un_xmax_min, un); un_xmax_max = max(un_xmax_max, un)
        end
    end
    if y_min_io
        @inbounds for k in 3:grid.mz-2, i in 3:grid.mx-2
            un = -buffers.v_face_y[i, 3, k]
            area = grid.dx * grid.dz[k]
            flow_ymin += un * area
            if un >= 0.0; flow_ymin_out += un * area; else; flow_ymin_in += -un * area; end
            un_ymin_min = min(un_ymin_min, un); un_ymin_max = max(un_ymin_max, un)
        end
    end
    if y_max_io
        @inbounds for k in 3:grid.mz-2, i in 3:grid.mx-2
            un = buffers.v_face_y[i, grid.my-1, k]
            area = grid.dx * grid.dz[k]
            flow_ymax += un * area
            if un >= 0.0; flow_ymax_out += un * area; else; flow_ymax_in += -un * area; end
            un_ymax_min = min(un_ymax_min, un); un_ymax_max = max(un_ymax_max, un)
        end
    end
    if z_min_io
        @inbounds for j in 3:grid.my-2, i in 3:grid.mx-2
            un = -buffers.w_face_z[i, j, 3]
            area = grid.dx * grid.dy
            flow_zmin += un * area
            if un >= 0.0; flow_zmin_out += un * area; else; flow_zmin_in += -un * area; end
            un_zmin_min = min(un_zmin_min, un); un_zmin_max = max(un_zmin_max, un)
        end
    end
    if z_max_io
        @inbounds for j in 3:grid.my-2, i in 3:grid.mx-2
            un = buffers.w_face_z[i, j, grid.mz-1]
            area = grid.dx * grid.dy
            flow_zmax += un * area
            if un >= 0.0; flow_zmax_out += un * area; else; flow_zmax_in += -un * area; end
            un_zmax_min = min(un_zmax_min, un); un_zmax_max = max(un_zmax_max, un)
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

function log_stability_violation(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    dt_fixed::Float64,
    metrics::NamedTuple, # from calculate_stability_metrics
    poisson_params::PoissonConfig,
    debug::Bool
)
    if debug
        cfl = metrics.cfl
        cfl_loc = metrics.cfl_loc
        u_max = metrics.u_max
        u_max_loc = metrics.u_max_loc
        cfl_i, cfl_j, cfl_k = cfl_loc
        u_i, u_j, u_k = u_max_loc
        
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
            if m0 == 0.0; continue; end
            
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
        
        # Poisson residuals
        update_boundary_mask!(buffers.mask, grid, bc_set, 0.0) # inflow/outflow/opening set to 0
        alpha = poisson_params.mach2 / (dt_fixed * dt_fixed)
        res0 = PressureSolver.compute_residual_sor(buffers.p_prev, buffers.rhs, buffers.mask, grid, poisson_params.omega, alpha)
        res1 = PressureSolver.compute_residual_sor(buffers.p, buffers.rhs, buffers.mask, grid, poisson_params.omega, alpha)
        update_boundary_mask!(buffers.mask, grid, bc_set, 1.0)
        
        println("Stability violation: CFL=$(metrics.cfl), D=$(metrics.diff_num)")
        println("  CFL max at (i,j,k)=($(cfl_i),$(cfl_j),$(cfl_k)) mask=$(m_c) u=$(u_c) v=$(v_c) w=$(w_c)")
        println("  Umax at (i,j,k)=($(u_i),$(u_j),$(u_k)) |U|=$(u_max) u=$(u_m) v=$(v_m) w=$(w_m)")
        println("  U* max at (i,j,k)=($(us_i),$(us_j),$(us_k)) |U*|=$(ustar_max)")
        println("  RHS max at (i,j,k)=($(rh_i),$(rh_j),$(rh_k)) |rhs|=$(rhs_max)")
        println("  |p| max at (i,j,k)=($(p_i),$(p_j),$(p_k)) |p|=$(p_max)")
        println("  |div(u*)| max at (i,j,k)=($(du_i),$(du_j),$(du_k)) |div|=$(divu_max)")
        println("  Poisson residuals: r0=$(res0) r1=$(res1)")
    else
        println("Stability violation: CFL=$(metrics.cfl), D=$(metrics.diff_num)")
    end
end

end # module Monitor
