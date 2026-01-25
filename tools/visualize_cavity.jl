#!/usr/bin/env julia
"""
Cavity Flow Visualization Utility

Reads SPH velocity field files and plots:
- w-velocity profile along x=0.5 (vertical centerline)
- u-velocity profile along z=0.5 (horizontal centerline)

Usage:
    julia tools/visualize_cavity.jl <sph_file>
    julia tools/visualize_cavity.jl <sph_file> <step_start> <step_end>
    julia tools/visualize_cavity.jl --config visualize.json

Example:
    julia tools/visualize_cavity.jl verification/cavity/output/vel_0001000.sph
    julia tools/visualize_cavity.jl verification/cavity/output/vel 1000 2000
"""

using CairoMakie
using JSON3
using Printf

include("sph_reader.jl")

struct SliceConfig
    plane::Symbol
    plane_index::Int
    variables::Vector{Symbol}
    output_format::Symbol
    output_dir::String
    vector_enabled::Bool
    vector_skip::Int
    text_output::Bool
    vector_scale::Float64
end

function infer_prs_path(vel_path::String)
    base = basename(vel_path)
    prs_base = startswith(base, "vel") ? replace(base, "vel" => "prs") : replace(base, "vel" => "prs", count=1)
    return joinpath(dirname(vel_path), prs_base)
end

function parse_slice_config(viz, options, base_dir::String)
    plane = Symbol(lowercase(String(get(viz, :plane, "xy"))))
    plane_index = Int(get(viz, :plane_index, 1))
    variables = [Symbol(lowercase(String(v))) for v in get(viz, :variables, ["velocity", "pressure"])]
    output_format = Symbol(lowercase(String(get(viz, :output_format, "png"))))
    output_dir = String(get(viz, :output_dir, "viz"))
    if !isabspath(output_dir)
        output_dir = joinpath(base_dir, output_dir)
    end
    vector_enabled = Bool(get(viz, :vector_enabled, false))
    vector_skip = Int(get(viz, :vector_skip, 1))
    text_output = Bool(get(viz, :text_output, false))
    vector_scale = isnothing(options) ? 0.05 : Float64(get(options, :vecscale, 0.05))
    return SliceConfig(plane, plane_index, variables, output_format, output_dir, vector_enabled, vector_skip, text_output, vector_scale)
end

function write_slice_text(filepath::String, vd, pd, config::SliceConfig)
    mkpath(dirname(filepath))
    x = [vd.x0 + (i - 0.5) * vd.dx for i in 1:vd.nx]
    y = [vd.y0 + (j - 0.5) * vd.dy for j in 1:vd.ny]
    z = [vd.z0 + (k - 0.5) * vd.dz for k in 1:vd.nz]

    u = vd.u
    v = vd.v
    w = vd.w
    p = isnothing(pd) ? nothing : pd.p

    open(filepath, "w") do io
        println(io, "x y z u v w p")
        if config.plane == :xy
            k = clamp(config.plane_index, 1, vd.nz)
            for j in 1:vd.ny, i in 1:vd.nx
                pval = isnothing(p) ? 0.0 : p[i, j, k]
                @printf(io, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                        x[i], y[j], z[k], u[i, j, k], v[i, j, k], w[i, j, k], pval)
            end
        elseif config.plane == :xz
            j = clamp(config.plane_index, 1, vd.ny)
            for k in 1:vd.nz, i in 1:vd.nx
                pval = isnothing(p) ? 0.0 : p[i, j, k]
                @printf(io, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                        x[i], y[j], z[k], u[i, j, k], v[i, j, k], w[i, j, k], pval)
            end
        elseif config.plane == :yz
            i = clamp(config.plane_index, 1, vd.nx)
            for k in 1:vd.nz, j in 1:vd.ny
                pval = isnothing(p) ? 0.0 : p[i, j, k]
                @printf(io, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                        x[i], y[j], z[k], u[i, j, k], v[i, j, k], w[i, j, k], pval)
            end
        end
    end
end

function render_slice(vd, pd, config::SliceConfig)
    if config.output_format == :none && !config.text_output
        return
    end

    x = [vd.x0 + (i - 0.5) * vd.dx for i in 1:vd.nx]
    y = [vd.y0 + (j - 0.5) * vd.dy for j in 1:vd.ny]
    z = [vd.z0 + (k - 0.5) * vd.dz for k in 1:vd.nz]

    u = vd.u
    v = vd.v
    w = vd.w
    p = isnothing(pd) ? nothing : pd.p

    plot_vars = [v for v in config.variables if v == :velocity || v == :pressure]
    if isempty(plot_vars)
        return
    end

    mkpath(config.output_dir)
    filename = joinpath(config.output_dir, "viz_$(lpad(vd.step, 7, '0')).$(config.output_format)")
    text_path = joinpath(config.output_dir, "slice_$(lpad(vd.step, 7, '0')).txt")

    if config.text_output
        write_slice_text(text_path, vd, pd, config)
    end
    if config.output_format == :none
        return
    end

    nplots = length(plot_vars)
    fig = Figure(size = (500 * nplots, 500))

    if config.plane == :xy
        k = clamp(config.plane_index, 1, vd.nz)
        X = x
        Y = y
        Z_val = z[k]

        data_mag = sqrt.(u[:, :, k].^2 .+ v[:, :, k].^2 .+ w[:, :, k].^2)
        data_u = u[:, :, k]
        data_v = v[:, :, k]
        data_p = isnothing(p) ? nothing : p[:, :, k]

        for (col, vname) in enumerate(plot_vars)
            if vname == :velocity
                ax = Axis(fig[1, 2*col-1], title = "Velocity Mag. Z=$(round(Z_val, digits=3))",
                          xlabel="X [m]", ylabel="Y [m]")
                hm = heatmap!(ax, X, Y, data_mag, colormap = :viridis)
                Colorbar(fig[1, 2*col], hm, label = "Vel [m/s]")
                if config.vector_enabled
                    skip = config.vector_skip
                    arrows2d!(ax, X[1:skip:end], Y[1:skip:end],
                              data_u[1:skip:end, 1:skip:end], data_v[1:skip:end, 1:skip:end],
                              tipwidth=7.5, tiplength=7.5, lengthscale=config.vector_scale, color=:white)
                end
            elseif vname == :pressure && !isnothing(data_p)
                ax = Axis(fig[1, 2*col-1], title = "Pressure Z=$(round(Z_val, digits=3))",
                          xlabel="X [m]", ylabel="Y [m]")
                hm = heatmap!(ax, X, Y, data_p, colormap = :plasma)
                contour!(ax, X, Y, data_p, color = :black, linewidth = 0.5, alpha = 0.5)
                Colorbar(fig[1, 2*col], hm, label = "P")
            end
        end
    elseif config.plane == :xz
        j = clamp(config.plane_index, 1, vd.ny)
        X = x
        Z = z
        Y_val = y[j]

        data_mag = sqrt.(u[:, j, :].^2 .+ v[:, j, :].^2 .+ w[:, j, :].^2)
        data_u = u[:, j, :]
        data_w = w[:, j, :]
        data_p = isnothing(p) ? nothing : p[:, j, :]

        for (col, vname) in enumerate(plot_vars)
            if vname == :velocity
                ax = Axis(fig[1, 2*col-1], title = "Velocity Mag. Y=$(round(Y_val, digits=3))",
                          xlabel="X [m]", ylabel="Z [m]")
                hm = heatmap!(ax, X, Z, data_mag, colormap = :viridis)
                Colorbar(fig[1, 2*col], hm, label = "Vel [m/s]")
                if config.vector_enabled
                    skip = config.vector_skip
                    arrows2d!(ax, X[1:skip:end], Z[1:skip:end],
                              data_u[1:skip:end, 1:skip:end], data_w[1:skip:end, 1:skip:end],
                              tipwidth=7.5, tiplength=7.5, lengthscale=config.vector_scale, color=:white)
                end
            elseif vname == :pressure && !isnothing(data_p)
                ax = Axis(fig[1, 2*col-1], title = "Pressure Y=$(round(Y_val, digits=3))",
                          xlabel="X [m]", ylabel="Z [m]")
                hm = heatmap!(ax, X, Z, data_p, colormap = :plasma)
                contour!(ax, X, Z, data_p, color = :black, linewidth = 0.5, alpha = 0.5)
                Colorbar(fig[1, 2*col], hm, label = "P")
            end
        end
    elseif config.plane == :yz
        i = clamp(config.plane_index, 1, vd.nx)
        Y = y
        Z = z
        X_val = x[i]

        data_mag = sqrt.(u[i, :, :].^2 .+ v[i, :, :].^2 .+ w[i, :, :].^2)
        data_v = v[i, :, :]
        data_w = w[i, :, :]
        data_p = isnothing(p) ? nothing : p[i, :, :]

        for (col, vname) in enumerate(plot_vars)
            if vname == :velocity
                ax = Axis(fig[1, 2*col-1], title = "Velocity Mag. X=$(round(X_val, digits=3))",
                          xlabel="Y [m]", ylabel="Z [m]")
                hm = heatmap!(ax, Y, Z, data_mag, colormap = :viridis)
                Colorbar(fig[1, 2*col], hm, label = "Vel [m/s]")
                if config.vector_enabled
                    skip = config.vector_skip
                    arrows2d!(ax, Y[1:skip:end], Z[1:skip:end],
                              data_v[1:skip:end, 1:skip:end], data_w[1:skip:end, 1:skip:end],
                              tipwidth=7.5, tiplength=7.5, lengthscale=config.vector_scale, color=:white)
                end
            elseif vname == :pressure && !isnothing(data_p)
                ax = Axis(fig[1, 2*col-1], title = "Pressure X=$(round(X_val, digits=3))",
                          xlabel="Y [m]", ylabel="Z [m]")
                hm = heatmap!(ax, Y, Z, data_p, colormap = :plasma)
                contour!(ax, Y, Z, data_p, color = :black, linewidth = 0.5, alpha = 0.5)
                Colorbar(fig[1, 2*col], hm, label = "P")
            end
        end
    end

    save(filename, fig)
    println("Saved: $filename")
end

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

function render_profile_from_data(vd, vel_path::String; output_dir::Union{Nothing,String}=nothing)
    profiles = extract_centerline_profiles(vd)
    println("  Y plane: j=$(profiles.j_mid) (y=$(round(profiles.y_val, digits=3)))")
    println("  X line:  i=$(profiles.i_mid) (x=$(round(profiles.x_val, digits=3)))")
    println("  Z line:  k=$(profiles.k_mid) (z=$(round(profiles.z_val, digits=3)))")

    basename_noext = replace(basename(vel_path), r"\.sph$" => "")
    output_dir = isnothing(output_dir) ? dirname(vel_path) : output_dir
    mkpath(output_dir)
    output_path = joinpath(output_dir, "$(basename_noext)_profile.png")

    plot_profiles(vd, profiles, output_path=output_path)
end

function run_slice_from_file(vel_path::String, config::SliceConfig; prs_path::Union{Nothing,String}=nothing)
    vd = read_sph_vector(vel_path)
    need_prs = any(v -> v == :pressure, config.variables)
    pd = nothing
    if need_prs
        p_path = isnothing(prs_path) ? infer_prs_path(vel_path) : prs_path
        if !isfile(p_path)
            error("Pressure file not found: $(p_path)")
        end
        pd = read_sph_scalar(p_path)
    end
    render_slice(vd, pd, config)
    render_profile_from_data(vd, vel_path; output_dir=config.output_dir)
end

function run_slice_from_range(prefix::String, step_start::Int, step_end::Int, config::SliceConfig)
    for step in step_start:step_end
        vel_path = "$(prefix)_$(lpad(step, 7, '0')).sph"
        if isfile(vel_path)
            run_slice_from_file(vel_path, config)
        end
    end
end

function run_from_config(path::String)
    if !isfile(path)
        error("Config file not found: $(path)")
    end
    data = JSON3.read(read(path, String))
    config_dir = dirname(path)
    mode = lowercase(String(get(data, :mode, "profile")))
    input = get(data, :input, nothing)
    if isnothing(input)
        error("Config must include input")
    end
    resolve_path(p::String) = isabspath(p) ? p : joinpath(config_dir, p)
    vel_path = haskey(input, :vel) ? resolve_path(String(input.vel)) : nothing
    sph_path = haskey(input, :sph) ? resolve_path(String(input.sph)) : nothing
    prefix_path = haskey(input, :prefix) ? resolve_path(String(input.prefix)) : nothing
    base_dir = !isnothing(vel_path) ? dirname(vel_path) :
               !isnothing(sph_path) ? dirname(sph_path) :
               !isnothing(prefix_path) ? dirname(prefix_path) : "."

    if mode == "slice"
        viz = haskey(data, :Visualization) ? data.Visualization : get(data, :visualization, nothing)
        if isnothing(viz)
            error("Slice mode requires visualization settings")
        end
        options = get(data, :options, nothing)
        config = parse_slice_config(viz, options, base_dir)

        if !isnothing(vel_path)
            prs_path = haskey(input, :prs) ? resolve_path(String(input.prs)) : nothing
            run_slice_from_file(vel_path, config; prs_path=prs_path)
            return
        end
        if !isnothing(prefix_path) && haskey(input, :step_start) && haskey(input, :step_end)
            run_slice_from_range(prefix_path, Int(input.step_start), Int(input.step_end), config)
            return
        end
        error("Slice mode input must include vel or (prefix, step_start, step_end)")
    else
        if !isnothing(sph_path)
            process_file(sph_path)
            return
        end
        if !isnothing(vel_path)
            process_file(vel_path)
            return
        end
        if !isnothing(prefix_path) && haskey(input, :step_start) && haskey(input, :step_end)
            process_range(prefix_path, Int(input.step_start), Int(input.step_end))
            return
        end
        error("Profile mode input must include sph/vel or (prefix, step_start, step_end)")
    end
end

# Main entry point
function main()
    if length(ARGS) < 1
        println("Usage:")
        println("  julia visualize_cavity.jl <sph_file>")
        println("  julia visualize_cavity.jl <prefix> <step_start> <step_end>")
        println("  julia visualize_cavity.jl --config visualize.json")
        return
    end

    args = collect(ARGS)
    if length(args) >= 2 && args[1] == "--config"
        run_from_config(args[2])
        return
    elseif length(args) >= 1 && startswith(args[1], "--config=")
        run_from_config(split(args[1], "=", limit=2)[2])
        return
    end

    if length(args) == 1
        # Single file
        process_file(args[1])
    elseif length(args) >= 3
        # Range of files
        prefix = args[1]
        step_start = parse(Int, args[2])
        step_end = parse(Int, args[3])
        process_range(prefix, step_start, step_end)
    else
        println("Error: Invalid arguments")
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
