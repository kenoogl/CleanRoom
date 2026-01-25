#!/usr/bin/env julia
"""
Cavity Flow Visualization Utility

Reads SPH velocity field files and plots:
- w-velocity profile along x=0.5 (vertical centerline)
- u-velocity profile along z=0.5 (horizontal centerline)

Usage:
    julia tools/visualize_cavity.jl <sph_file>
    julia tools/visualize_cavity.jl <sph_file> <step_start> <step_end>

Example:
    julia tools/visualize_cavity.jl verification/cavity/output/vel_0001000.sph
    julia tools/visualize_cavity.jl verification/cavity/output/vel 1000 2000
"""

using CairoMakie

include("sph_reader.jl")

"""
Extract centerline profiles at Y=0.5 plane for Ghia-style plots
- u along vertical centerline (x=0.5): u(z)
- w along horizontal centerline (z=0.5): w(x)
"""
function extract_centerline_profiles(data)
    nx, ny, nz = data.nx, data.ny, data.nz
    x0, y0, z0 = data.x0, data.y0, data.z0
    dx, dy, dz = data.dx, data.dy, data.dz
    
    # Coordinate arrays (cell centers)
    # SPH origin is at left edge of first cell
    # Cell center = origin + (i - 0.5) * dx
    x = [x0 + (i - 0.5) * dx for i in 1:nx]
    y = [y0 + (j - 0.5) * dy for j in 1:ny]
    z = [z0 + (k - 0.5) * dz for k in 1:nz]
    
    # Find Y = 0.5 plane index
    y_target = 0.5
    j_mid = argmin(abs.(y .- y_target))
    
    # Find X = 0.5 index (for u profile along z - vertical centerline)
    x_target = 0.5
    i_mid = argmin(abs.(x .- x_target))
    
    # Find Z = 0.5 index (for w profile along x - horizontal centerline)
    z_target = 0.5
    k_mid = argmin(abs.(z .- z_target))
    
    # Extract profiles (Ghia-style)
    # u along vertical centerline (x=0.5): u(z) at i_mid, j_mid
    u_centerline = data.u[i_mid, j_mid, :]
    
    # w along horizontal centerline (z=0.5): w(x) at j_mid, k_mid
    w_centerline = data.w[:, j_mid, k_mid]
    
    return (z_coords=z, u_centerline=u_centerline,
            x_coords=x, w_centerline=w_centerline,
            j_mid=j_mid, i_mid=i_mid, k_mid=k_mid,
            y_val=y[j_mid], x_val=x[i_mid], z_val=z[k_mid])
end

"""
Plot cavity flow centerline profiles in Ghia et al. paper format
- Left panel: u velocity vs z (vertical centerline at x=0.5)
- Right panel: w velocity vs x (horizontal centerline at z=0.5)
"""
function plot_profiles(data, profiles; output_path=nothing, title_suffix="")
    # Two-panel layout
    fig = Figure(size = (800, 400))
    
    supertitle = "Cavity Flow (Step $(data.step), t=$(round(data.time, digits=2)))$(title_suffix)"
    Label(fig[0, :], supertitle, fontsize = 20)
    
    # Left panel: u vs z (center u-axis at zero)
    max_u = maximum(abs.(profiles.u_centerline))
    u_lim = max(max_u * 1.05, 1.0e-6)
    ax1 = Axis(fig[1, 1], xlabel = "u", ylabel = "z", 
               title = "u at x=$(round(profiles.x_val, digits=2))",
               limits = ((-u_lim, u_lim), (0, 1)),
               aspect = AxisAspect(1),
               xticklabelpad = 5, yticklabelpad = 5,
               xlabelpadding = 10, ylabelpadding = 10)
    
    scatterlines!(ax1, profiles.u_centerline, profiles.z_coords, 
                  label = "Present", color = :blue, markersize = 8)
    axislegend(ax1, position = :rb)
    
    # Right panel: w vs x (center w-axis at zero)
    max_w = maximum(abs.(profiles.w_centerline))
    w_lim = max(max_w * 1.05, 1.0e-6)
    ax2 = Axis(fig[1, 2], xlabel = "x", ylabel = "w", 
               title = "w at z=$(round(profiles.z_val, digits=2))",
               limits = ((0, 1), (-w_lim, w_lim)),
               aspect = AxisAspect(1),
               xticklabelpad = 5, yticklabelpad = 5,
               xlabelpadding = 10, ylabelpadding = 10)
    
    scatterlines!(ax2, profiles.x_coords, profiles.w_centerline, 
                  label = "Present", color = :blue, markersize = 8)
    axislegend(ax2, position = :rt)
    
    if !isnothing(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end
    
    return fig
end

"""
Process single SPH file
"""
function process_file(filepath::String)
    println("Reading: $filepath")
    data = read_sph_vector(filepath)
    println("  Grid: $(data.nx) x $(data.ny) x $(data.nz)")
    println("  Step: $(data.step), Time: $(data.time)")
    
    profiles = extract_centerline_profiles(data)
    println("  Y plane: j=$(profiles.j_mid) (y=$(round(profiles.y_val, digits=3)))")
    println("  X line:  i=$(profiles.i_mid) (x=$(round(profiles.x_val, digits=3)))")
    println("  Z line:  k=$(profiles.k_mid) (z=$(round(profiles.z_val, digits=3)))")
    
    # Generate output filename
    basename_noext = replace(basename(filepath), r"\.sph$" => "")
    output_dir = dirname(filepath)
    output_path = joinpath(output_dir, "$(basename_noext)_profile.png")
    
    plot_profiles(data, profiles, output_path=output_path)
    
    return (data=data, profiles=profiles)
end

"""
Process range of SPH files
"""
function process_range(prefix::String, step_start::Int, step_end::Int)
    results = []
    
    for step in step_start:step_end
        filepath = "$(prefix)_$(lpad(step, 7, '0')).sph"
        if isfile(filepath)
            result = process_file(filepath)
            push!(results, result)
        end
    end
    
    return results
end

# Main entry point
function main()
    if length(ARGS) < 1
        println("Usage:")
        println("  julia visualize_cavity.jl <sph_file>")
        println("  julia visualize_cavity.jl <prefix> <step_start> <step_end>")
        return
    end
    
    if length(ARGS) == 1
        # Single file
        process_file(ARGS[1])
    elseif length(ARGS) >= 3
        # Range of files
        prefix = ARGS[1]
        step_start = parse(Int, ARGS[2])
        step_end = parse(Int, ARGS[3])
        process_range(prefix, step_start, step_end)
    else
        println("Error: Invalid arguments")
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
