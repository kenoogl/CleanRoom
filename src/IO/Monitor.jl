module Monitor

using Printf
using ..Common
using ..Fields
using ..Grid

export MonitorData, MonitorConfig, init_monitor, log_step!, check_divergence
export compute_u_max, compute_cfl, compute_divergence_max

struct MonitorData
    step::Int
    time::Float64
    u_max::Float64
    div_max::Float64
    dU::Float64  # Velocity change for steady-state check: sqrt(Σ((u^n+1-u^n)² + ...))
    pressure_itr::Int
    pressure_residual::Float64
end

struct MonitorConfig
    console_interval::Int
    history_interval::Int
    div_threshold::Float64
    step_digits::Int
end

function init_monitor(config::MonitorConfig, history_path::String)
    open(history_path, "w") do io
        # Header alignment to match data
        # Step: width=config.step_digits (left aligned)
        # Time: %.6e -> ~13 chars + 1 space -> 14
        # dt: %.5e -> ~12 chars + 1 space -> 13
        # Umax: %.4e -> ~11 chars + 1 space -> 12
        # CFL: %.2f -> ~6 chars + 1 space -> 7
        # divMax: %.4e -> ~11 chars + 1 space -> 12
        # ItrP: %4d -> 4 chars + 1 space -> 5
        # ResP: %.5e -> ~12 chars + 1 space -> 13
        
        # We construct header string manually or just use broad enough spacing
        # Let's use specific spacing
        # Ensure step width accommodates "step" (4 chars)
        step_width = max(config.step_digits, 4)
        s_step = rpad("step", step_width)
        
        # Using approximated spacing based on the printf below
        # "%-Nd %14.6e %13.5e %12.4e %7.2f %12.4e %5d %13.5e"
        
        print(io, s_step, " ")
        # Right aligned headers to match numeric data
        @printf(io, "%14s %12s %12s %12s %5s %13s\n",
            "time", "Umax", "divMax", "dU", "ItrP", "ResP")
    end
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

function compute_divergence_max(buffers::CFDBuffers, grid::GridData)
    max_div = 0.0
    dx, dy = grid.dx, grid.dy
    @inbounds for k in 3:grid.mz-2
        dz = grid.dz[k]
        for j in 3:grid.my-2, i in 3:grid.mx-2
            # Central diff for div u
             div = (buffers.u[i+1, j, k] - buffers.u[i-1, j, k])/(2*dx) +
                   (buffers.v[i, j+1, k] - buffers.v[i, j-1, k])/(2*dy) +
                   (buffers.w[i, j, k+1] - buffers.w[i, j, k-1])/(2*dz) 
             if abs(div) > max_div; max_div = abs(div); end
        end
    end
    return max_div
end

function log_step!(
    data::MonitorData,
    config::MonitorConfig,
    history_io::IO
)
    if data.step % config.console_interval == 0 || data.step == 1
        @printf("Step %d Time %.4e Umax %.3e Div %.3e dU %.3e ItrP %d Pres %.3e\n",
            data.step, data.time, data.u_max, data.div_max, data.dU, data.pressure_itr, data.pressure_residual)
    end
    
    if data.step % config.history_interval == 0 || data.step == 1
        # Step: Left aligned, padded to max digits
        step_width = max(config.step_digits, 4)
        s_step = rpad(string(data.step), step_width)
        print(history_io, s_step, " ")
        
        # Alignment:
        # time:   .6e (e.g. 1.234567e+00) -> 12 chars. Allow 14.
        # Umax:   .4e (e.g. 1.2345e+00) -> 10 chars. Allow 12.
        # divMax: .4e -> 10 chars. Allow 12.
        # ItrP:   %4d -> 4 chars. Allow 5.
        # ResP:   .5e -> 11 chars. Allow 13.
        
        @printf(history_io, "%14.6e %12.4e %12.4e %12.4e %5d %13.5e\n",
            data.time, data.u_max, data.div_max, data.dU, data.pressure_itr, data.pressure_residual)
        flush(history_io)
    end
end

function check_divergence(div_max::Float64, threshold::Float64)::Bool
    return div_max > threshold || isnan(div_max) || isinf(div_max)
end

end # module Monitor
