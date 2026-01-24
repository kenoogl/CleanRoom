module CleanroomVisualizationExt

using CairoMakie
using CleanroomSolver
using CleanroomSolver.Visualization
using CleanroomSolver.Fields
using CleanroomSolver.Grid
using CleanroomSolver.Common

# Helper to reconstruct coords
function get_coords(grid::GridData, dim_params::DimensionParams)
    # Dimensional coords
    x = grid.x .* dim_params.L0
    y = grid.y .* dim_params.L0
    z_center = grid.z_center .* dim_params.L0
    return x, y, z_center
end

@inline function clamp_internal_index(idx::Int, max_internal::Int)::Int
    return clamp(idx, 1, max_internal) + 2
end

# Extend render_slice
function Visualization.render_slice(
    buffers::CFDBuffers,
    grid::GridData,
    config::VizConfig,
    step::Int,
    dim_params::DimensionParams
)
    if config.text_output
        Visualization.write_slice_text("$(config.output_dir)/slice_$(step).txt", buffers, grid, config, dim_params)
    end

    if config.output_format == :none
        return
    end

    plot_vars = [v for v in config.variables if v == :velocity || v == :pressure]
    if isempty(plot_vars)
        return
    end

    mkpath(config.output_dir)
    filename = "$(config.output_dir)/viz_$(lpad(step, 7, '0')).$(config.output_format)"
    
    # 1. Extract Slice Data
    x, y, z = get_coords(grid, dim_params)
    U0 = dim_params.U0
    
    # Slice Indexing
    idx = config.plane_index

    # Prepare Data for Plotting
    # Needs to be a Matrix for Heatmap
    # X-axis, Y-axis, Values
    var_u = buffers.u .* U0
    var_v = buffers.v .* U0
    var_w = buffers.w .* U0
    var_p = buffers.p .* (U0^2)
    
    # Magnitude
    var_mag = sqrt.(var_u.^2 .+ var_v.^2 .+ var_w.^2)
    
    nplots = length(plot_vars)
    fig = Figure(size = (500 * nplots, 500))
    
    # Plot Logic based on plane
    if config.plane == :xy
        # XY plane at k = idx (internal index)
        k = clamp_internal_index(idx, grid.mz - 4)
        
        # Slicing (excluding ghosts if we want strictly domain, but buffers includes ghosts)
        # grid.x is full including ghosts.
        # We should probably plot valid domain.
        # Indices: 3 : mx-2
        valid_x = 3:grid.mx-2
        valid_y = 3:grid.my-2
        
        X = x[valid_x]
        Y = y[valid_y]
        Z_val = z[k]
        
        # Extract 2D array
        data_mag = var_mag[valid_x, valid_y, k]
        data_u   = var_u[valid_x, valid_y, k]
        data_v   = var_v[valid_x, valid_y, k]
        data_p   = var_p[valid_x, valid_y, k]
        
        for (col, vname) in enumerate(plot_vars)
            if vname == :velocity
                ax = Axis(fig[1, 2*col-1], title = "Velocity Mag. Z=$(round(Z_val, digits=3))m", xlabel="X [m]", ylabel="Y [m]")
                hm = heatmap!(ax, X, Y, data_mag, colormap = :viridis)
                Colorbar(fig[1, 2*col], hm, label = "Vel [m/s]")
                
                if config.vector_enabled
                    skip = config.vector_skip
                    arrows2d!(ax, X[1:skip:end], Y[1:skip:end], data_u[1:skip:end, 1:skip:end], data_v[1:skip:end, 1:skip:end],
                              tipwidth=7.5, tiplength=7.5, lengthscale=0.05 / U0, color=:white)
                end
            elseif vname == :pressure
                ax = Axis(fig[1, 2*col-1], title = "Pressure Z=$(round(Z_val, digits=3))m", xlabel="X [m]", ylabel="Y [m]")
                hm = heatmap!(ax, X, Y, data_p, colormap = :plasma)
                contour!(ax, X, Y, data_p, color = :black, linewidth = 0.5, alpha = 0.5)
                Colorbar(fig[1, 2*col], hm, label = "P [Pa]")
            end
        end
        
    elseif config.plane == :xz
        # XZ plane at j = idx
        j = clamp_internal_index(idx, grid.my - 4)
        
        valid_x = 3:grid.mx-2
        valid_z = 3:grid.mz-2
        
        X = x[valid_x]
        Z = z[valid_z]
        Y_val = y[j]
        
        data_mag = var_mag[valid_x, j, valid_z]
        data_u   = var_u[valid_x, j, valid_z]
        data_w   = var_w[valid_x, j, valid_z]
        data_p   = var_p[valid_x, j, valid_z]
        
        for (col, vname) in enumerate(plot_vars)
            if vname == :velocity
                ax = Axis(fig[1, 2*col-1], title = "Velocity Mag. Y=$(round(Y_val, digits=3))m", xlabel="X [m]", ylabel="Z [m]")
                hm = heatmap!(ax, X, Z, data_mag, colormap = :viridis)
                Colorbar(fig[1, 2*col], hm, label = "Vel [m/s]")
                
                if config.vector_enabled
                    skip = config.vector_skip
                    arrows2d!(ax, X[1:skip:end], Z[1:skip:end], data_u[1:skip:end, 1:skip:end], data_w[1:skip:end, 1:skip:end],
                              tipwidth=7.5, tiplength=7.5, lengthscale=0.05 / U0, color=:white)
                end
            elseif vname == :pressure
                ax = Axis(fig[1, 2*col-1], title = "Pressure Y=$(round(Y_val, digits=3))m", xlabel="X [m]", ylabel="Z [m]")
                hm = heatmap!(ax, X, Z, data_p, colormap = :plasma)
                contour!(ax, X, Z, data_p, color = :black, linewidth = 0.5, alpha = 0.5)
                Colorbar(fig[1, 2*col], hm, label = "P [Pa]")
            end
        end
        
    elseif config.plane == :yz
         # YZ plane at i = idx
        i = clamp_internal_index(idx, grid.mx - 4)
        
        valid_y = 3:grid.my-2
        valid_z = 3:grid.mz-2
        
        Y = y[valid_y]
        Z = z[valid_z]
        X_val = x[i]
        
        data_mag = var_mag[i, valid_y, valid_z]
        data_v   = var_v[i, valid_y, valid_z]
        data_w   = var_w[i, valid_y, valid_z]
        data_p   = var_p[i, valid_y, valid_z]
        
        for (col, vname) in enumerate(plot_vars)
            if vname == :velocity
                ax = Axis(fig[1, 2*col-1], title = "Velocity Mag. X=$(round(X_val, digits=3))m", xlabel="Y [m]", ylabel="Z [m]")
                hm = heatmap!(ax, Y, Z, data_mag, colormap = :viridis)
                Colorbar(fig[1, 2*col], hm, label = "Vel [m/s]")
                
                if config.vector_enabled
                    skip = config.vector_skip
                    arrows2d!(ax, Y[1:skip:end], Z[1:skip:end], data_v[1:skip:end, 1:skip:end], data_w[1:skip:end, 1:skip:end],
                              tipwidth=7.5, tiplength=7.5, lengthscale=0.05 / U0, color=:white)
                end
            elseif vname == :pressure
                ax = Axis(fig[1, 2*col-1], title = "Pressure X=$(round(X_val, digits=3))m", xlabel="Y [m]", ylabel="Z [m]")
                hm = heatmap!(ax, Y, Z, data_p, colormap = :plasma)
                contour!(ax, Y, Z, data_p, color = :black, linewidth = 0.5, alpha = 0.5)
                Colorbar(fig[1, 2*col], hm, label = "P [Pa]")
            end
        end
    end
    
    save(filename, fig)
    println("Saved visualization to $filename")
end

end # module
