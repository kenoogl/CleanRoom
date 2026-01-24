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
    u, v, w = buffers.u, buffers.v, buffers.w
    @inbounds for k in 1:grid.mz, j in 1:grid.my, i in 1:grid.mx
        mag = sqrt(u[i, j, k]^2 + v[i, j, k]^2 + w[i, j, k]^2)
        if mag > umax; umax = mag; end
    end
    return umax
end

function compute_cfl(buffers::CFDBuffers, grid::GridData, dt::Float64)
    cfl_max = 0.0
    u, v, w = buffers.u, buffers.v, buffers.w
    dx, dy = grid.dx, grid.dy
    @inbounds for k in 3:grid.mz-2
        dz = grid.dz[k]
        for j in 3:grid.my-2, i in 3:grid.mx-2
             val = abs(u[i, j, k])/dx + abs(v[i, j, k])/dy + abs(w[i, j, k])/dz
             if val > cfl_max; cfl_max = val; end
        end
    end
    return cfl_max * dt
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

end # module Monitor
