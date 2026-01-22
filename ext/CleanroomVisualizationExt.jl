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

    mkpath(config.output_dir)
    filename = "$(config.output_dir)/viz_$(lpad(step, 7, '0')).$(config.output_format)"
    
    # 1. Extract Slice Data
    x, y, z = get_coords(grid, dim_params)
    U0 = dim_params.U0
    
    # Slice Indexing
    idx = config.plane_index
    if idx < 1 
        idx = 1
    end
    
    # Prepare Data for Plotting
    # Needs to be a Matrix for Heatmap
    # X-axis, Y-axis, Values
    var_u = buffers.u .* U0
    var_v = buffers.v .* U0
    var_w = buffers.w .* U0
    var_p = buffers.p .* (U0^2 * 1.225) # Approx Pa if Rho=1.225? Or just dimensionless p * U0^2 if density is normalized
    # Let's assume just p * U0^2 for now or "Pressure [Pa]" if rho included. 
    # Spec says "output_dimensional".
    
    # Magnitude
    var_mag = sqrt.(var_u.^2 .+ var_v.^2 .+ var_w.^2)
    
    title_str = "Step $step"
    
    fig = Figure(size = (1200, 500))
    
    # Plot Logic based on plane
    if config.plane == :xy
        # XY plane at k = idx
        if idx > grid.mz; idx = grid.mz; end
        
        # Slicing (excluding ghosts if we want strictly domain, but buffers includes ghosts)
        # grid.x is full including ghosts.
        # We should probably plot valid domain.
        # Indices: 3 : mx-2
        valid_x = 3:grid.mx-2
        valid_y = 3:grid.my-2
        
        X = x[valid_x]
        Y = y[valid_y]
        Z_val = z[idx]
        
        # Extract 2D array
        data_mag = var_mag[valid_x, valid_y, idx]
        data_u   = var_u[valid_x, valid_y, idx]
        data_v   = var_v[valid_x, valid_y, idx]
        
        ax = Axis(fig[1, 1], title = "Velocity Magnitude at Z=$(round(Z_val, digits=3)) m", xlabel="X [m]", ylabel="Y [m]", aspect=DataAspect())
        hm = heatmap!(ax, X, Y, data_mag, colormap = :viridis)
        Colorbar(fig[1, 2], hm, label = "Velocity [m/s]")
        
        if config.vector_enabled
             # Downsample for arrows
             skip = config.vector_skip
             arrows!(ax, X[1:skip:end], Y[1:skip:end], data_u[1:skip:end, 1:skip:end], data_v[1:skip:end, 1:skip:end], 
                     arrowsize=10, lengthscale=0.05 / U0, color=:white)
        end
        
    elseif config.plane == :xz
        # XZ plane at j = idx
        if idx > grid.my; idx = grid.my; end
        
        valid_x = 3:grid.mx-2
        valid_z = 3:grid.mz-2
        
        X = x[valid_x]
        Z = z[valid_z]
        Y_val = y[idx]
        
        data_mag = var_mag[valid_x, idx, valid_z]
        data_u   = var_u[valid_x, idx, valid_z]
        data_w   = var_w[valid_x, idx, valid_z]
        
        ax = Axis(fig[1, 1], title = "Velocity Magnitude at Y=$(round(Y_val, digits=3)) m", xlabel="X [m]", ylabel="Z [m]", aspect=DataAspect())
        hm = heatmap!(ax, X, Z, data_mag, colormap = :viridis)
        Colorbar(fig[1, 2], hm, label = "Velocity [m/s]")
        
        if config.vector_enabled
             skip = config.vector_skip
             arrows!(ax, X[1:skip:end], Z[1:skip:end], data_u[1:skip:end, 1:skip:end], data_w[1:skip:end, 1:skip:end], 
                     arrowsize=10, lengthscale=0.05 / U0, color=:white)
        end
        
    elseif config.plane == :yz
         # YZ plane at i = idx
        if idx > grid.mx; idx = grid.mx; end
        
        valid_y = 3:grid.my-2
        valid_z = 3:grid.mz-2
        
        Y = y[valid_y]
        Z = z[valid_z]
        X_val = x[idx]
        
        data_mag = var_mag[idx, valid_y, valid_z]
        data_v   = var_v[idx, valid_y, valid_z]
        data_w   = var_w[idx, valid_y, valid_z]
        
        ax = Axis(fig[1, 1], title = "Velocity Magnitude at X=$(round(X_val, digits=3)) m", xlabel="Y [m]", ylabel="Z [m]", aspect=DataAspect())
        hm = heatmap!(ax, Y, Z, data_mag, colormap = :viridis)
        Colorbar(fig[1, 2], hm, label = "Velocity [m/s]")
        
        if config.vector_enabled
             skip = config.vector_skip
             arrows!(ax, Y[1:skip:end], Z[1:skip:end], data_v[1:skip:end, 1:skip:end], data_w[1:skip:end, 1:skip:end], 
                     arrowsize=10, lengthscale=0.05 / U0, color=:white)
        end
    end
    
    save(filename, fig)
    println("Saved visualization to $filename")
end

end # module
