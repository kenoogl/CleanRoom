#!/usr/bin/env julia
using CairoMakie
using JSON3

include("sph_reader.jl")

"""
    plot_step_velocity(vd; output_path=nothing, vector_scale=0.05, vector_ref=nothing)

Plot velocity vectors for the backward facing step flow.
"""
function plot_step_velocity(vd; output_path=nothing, vector_scale=0.05, vector_ref=nothing)
    # Y-slice index (fixed as per original script, or could be dynamic)
    # Original used j=23. We can try to find the mid-plane if we know the domain,
    # but sticking to j=23 or middle index is safer if dimensions match.
    # Let's use the middle index if j=23 is not guaranteed, or keep 23 if it's specific.
    # The original script hardcoded j=23. Let's assume the grid hasn't changed drastically
    # or calculate the middle.
    j = div(vd.ny, 2) # Use approximate middle instead of hardcoded 23
    
    x = [vd.x0 + (i - 0.5) * vd.dx for i in 1:vd.nx]
    z = [vd.z0 + (k - 0.5) * vd.dz for k in 1:vd.nz]
    
    u = vd.u[:, j, :]
    w = vd.w[:, j, :]
    xrange = x[end] - x[1]
    zrange = z[end] - z[1]
    ratio = xrange / zrange
    display_ratio = ratio
    base_width = 1000
    title_height = 50
    plot_height = Int(round(base_width / display_ratio))
    fig = Figure(size = (base_width, plot_height + title_height), figure_padding = (10, 10, 40, 10))
    supertitle = "Step Flow Velocity (Step $(vd.step), t=$(round(vd.time, digits=2)))"
    Label(fig[1, :], supertitle, fontsize = 20)
    
    # Velocity Vectors (Arrows)
    zmax = z[1] + 1.2 * zrange
    ax1 = Axis(fig[2, 1], title = "Velocity Vectors", xlabel = "x [m]", ylabel = "z [m]",
               xlabelpadding = 12, xticklabelpad = 6,
               aspect = AxisAspect(display_ratio), limits = ((x[1], x[end]), (z[1], zmax)),
               xautolimitmargin = (0.0, 0.0), yautolimitmargin = (0.0, 0.0))
    
    # Downsample for visibility
    step_vec = 4
    arrows2d!(ax1, x[1:step_vec:end], z[1:step_vec:end],
              u[1:step_vec:end, 1:step_vec:end], w[1:step_vec:end, 1:step_vec:end],
              tipwidth = 7.5, tiplength = 7.5, lengthscale = vector_scale, color = :black)

    # Vector scale arrow (inside plot, above domain)
    speed = sqrt.(u.^2 .+ w.^2)
    max_speed = maximum(speed)
    ref_speed = isnothing(vector_ref) ? round(max_speed, sigdigits=3) : vector_ref
    if ref_speed > 0
        legend_len = ref_speed * vector_scale
        x0 = x[1] + 0.05 * xrange
        z0 = z[1] + 1.15 * zrange
        arrows2d!(ax1, [x0], [z0], [ref_speed], [0.0],
                  tipwidth = 7.5, tiplength = 7.5, lengthscale = vector_scale, color = :red)
        text!(ax1, "U = $(ref_speed)", position = (x0 + legend_len * 1.05, z0),
              align = (:left, :center), color = :red)
    end

    # Spacer column to match pressure plot width
    ax_spacer = Axis(fig[2, 2])
    hidedecorations!(ax_spacer)
    hidespines!(ax_spacer)

    rowsize!(fig.layout, 1, Fixed(title_height))
    rowsize!(fig.layout, 2, Fixed(plot_height))
    colsize!(fig.layout, 1, Relative(0.92))
    colsize!(fig.layout, 2, Relative(0.08))
    
    # Geometry Overlay
    # Step boundary (approximate based on standard BFS geometry)
    lines!(ax1, [-0.5, 0.0, 0.0], [0.3, 0.3, 0.0], color = :blue, linewidth = 2)
    # Solid fill
    poly!(ax1, Point2f[(-0.5, 0.0), (0.0, 0.0), (0.0, 0.3), (-0.5, 0.3)],
          color = (:gray, 0.3), strokewidth = 0)
    
    if !isnothing(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end
    
    return fig
end

"""
    plot_step_pressure(vd, pd; output_path=nothing)

Plot pressure distribution with contour overlay for the backward facing step flow.
"""
function plot_step_pressure(vd, pd; output_path=nothing)
    j = div(vd.ny, 2)

    x = [vd.x0 + (i - 0.5) * vd.dx for i in 1:vd.nx]
    z = [vd.z0 + (k - 0.5) * vd.dz for k in 1:vd.nz]

    p = pd.p[:, j, :]
    pmin, pmax = extrema(p)
    levels = range(pmin, pmax, length = 21) # 20 equal intervals

    xrange = x[end] - x[1]
    zrange = z[end] - z[1]
    ratio = xrange / zrange
    display_ratio = ratio
    base_width = 1000
    title_height = 50
    plot_height = Int(round(base_width / display_ratio))
    fig = Figure(size = (base_width, plot_height + title_height), figure_padding = (10, 10, 40, 10))
    supertitle = "Step Flow Pressure (Step $(vd.step), t=$(round(vd.time, digits=2)))"
    Label(fig[1, :], supertitle, fontsize = 20)

    ax = Axis(fig[2, 1], title = "Pressure (Heatmap + Contours)", xlabel = "x [m]", ylabel = "z [m]",
              xlabelpadding = 12, xticklabelpad = 6,
              aspect = AxisAspect(display_ratio), limits = ((x[1], x[end]), (z[1], z[end])),
              xautolimitmargin = (0.0, 0.0), yautolimitmargin = (0.0, 0.0))
    hm = heatmap!(ax, x, z, p, colormap = :plasma)
    contour!(ax, x, z, p, levels = levels, color = :black, linewidth = 0.8)
    Colorbar(fig[2, 2], hm, label = "P")

    colsize!(fig.layout, 1, Relative(0.92))
    colsize!(fig.layout, 2, Relative(0.08))
    rowsize!(fig.layout, 1, Fixed(title_height))
    rowsize!(fig.layout, 2, Fixed(plot_height))

    lines!(ax, [-0.5, 0.0, 0.0], [0.3, 0.3, 0.0], color = :blue, linewidth = 2)
    poly!(ax, Point2f[(-0.5, 0.0), (0.0, 0.0), (0.0, 0.3), (-0.5, 0.3)],
          color = (:gray, 0.3), strokewidth = 0)

    if !isnothing(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end

    return fig
end

"""
    process_file(vel_path; vector_scale=vector_scale, output_dir=nothing)

Process a single velocity file (and corresponding pressure file).
"""
function process_file(vel_path::String; vector_scale=0.05, vector_ref=nothing, output_dir::Union{Nothing,String}=nothing)
    println("Processing: $vel_path")
    
    # Infer pressure file path
    # Assume naming convention: vel_XXXXXXX.sph -> prs_XXXXXXX.sph
    dir = dirname(vel_path)
    out_dir = isnothing(output_dir) ? dir : output_dir
    mkpath(out_dir)
    base = basename(vel_path)
    
    # Replace "vel" with "prs" in the filename
    if startswith(base, "vel")
        prs_base = replace(base, "vel" => "prs")
    else
        # Fallback or error? Let's try to replace first occurrence of vel
        prs_base = replace(base, "vel" => "prs", count=1)
    end
    
    prs_path = joinpath(dir, prs_base)
    
    if !isfile(vel_path)
        println("Error: Velocity file not found: $vel_path")
        return
    end
    
    if !isfile(prs_path)
        println("Error: Pressure file not found: $prs_path")
        println("  Expected: $prs_path")
        return
    end
    
    vd = read_sph_vector(vel_path)
    pd = read_sph_scalar(prs_path)
    
    # Generate output path
    # Output: vel_XXXXXXX_plot.png
    out_base = replace(base, ".sph" => "_plot.png")
    out_path = joinpath(dir, out_base)
    
    step_match = match(r"vel_(\d+)\.sph$", base)
    if step_match !== nothing
        step_tag = step_match.captures[1]
        vel_out = joinpath(out_dir, "velocity_$(step_tag).png")
        prs_out = joinpath(out_dir, "pressure_$(step_tag).png")
    else
        vel_out = joinpath(out_dir, replace(base, ".sph" => "_velocity.png"))
        prs_out = joinpath(out_dir, replace(base, ".sph" => "_pressure.png"))
    end

    plot_step_velocity(vd, output_path=vel_out, vector_scale=vector_scale, vector_ref=vector_ref)
    plot_step_pressure(vd, pd, output_path=prs_out)
end

"""
    process_range(prefix, step_start, step_end)

Process a range of steps.
"""
function process_range(prefix::String, step_start::Int, step_end::Int; vector_scale=0.05, vector_ref=nothing, output_dir::Union{Nothing,String}=nothing)
    for step in step_start:step_end
        # Construct filename. Prefix is likely "verification/backward_step/output/vel"
        # We need to append "_XXXXXXX.sph"
        
        # Check if prefix already includes underscore or not
        # Standard usage: tools/visualize_step.jl verification/backward_step/output/vel 100 200
        
        fname = "$(prefix)_$(lpad(step, 7, '0')).sph"
        if isfile(fname)
            process_file(fname; vector_scale=vector_scale, vector_ref=vector_ref, output_dir=output_dir)
        end
    end
end

function parse_vector_scale(args::Vector{String})
    vector_scale = 0.05
    vector_ref = nothing
    filtered = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if startswith(arg, "--vecscale=")
            value = split(arg, "=", limit=2)[2]
            vector_scale = parse(Float64, value)
            i += 1
        elseif startswith(arg, "--vecref=")
            value = split(arg, "=", limit=2)[2]
            vector_ref = parse(Float64, value)
            i += 1
        elseif arg == "--vecscale"
            if i == length(args)
                error("Missing value for --vecscale")
            end
            vector_scale = parse(Float64, args[i + 1])
            i += 2
        elseif arg == "--vecref"
            if i == length(args)
                error("Missing value for --vecref")
            end
            vector_ref = parse(Float64, args[i + 1])
            i += 2
        else
            push!(filtered, arg)
            i += 1
        end
    end
    return vector_scale, vector_ref, filtered
end

function run_from_config(path::String)
    if !isfile(path)
        error("Config file not found: $(path)")
    end
    data = JSON3.read(read(path, String))
    input = get(data, :input, nothing)
    if isnothing(input)
        error("Config must include input")
    end
    config_dir = dirname(path)
    resolve_path(p::String) = isabspath(p) ? p : joinpath(config_dir, p)
    options = get(data, :options, nothing)
    vecscale = isnothing(options) ? 0.05 : Float64(get(options, :vecscale, 0.05))
    vecref = isnothing(options) ? nothing : (haskey(options, :vecref) ? Float64(options.vecref) : nothing)
    viz = haskey(data, :Visualization) ? data.Visualization : get(data, :visualization, nothing)
    outdir_opt = if !isnothing(options) && haskey(options, :output_dir)
        String(options.output_dir)
    elseif !isnothing(viz) && haskey(viz, :output_dir)
        String(viz.output_dir)
    else
        nothing
    end
    output_dir = isnothing(outdir_opt) ? nothing : resolve_path(outdir_opt)

    if haskey(input, :vel)
        process_file(resolve_path(String(input.vel)); vector_scale=vecscale, vector_ref=vecref, output_dir=output_dir)
        return
    end
    if haskey(input, :prefix) && haskey(input, :step_start) && haskey(input, :step_end)
        process_range(resolve_path(String(input.prefix)), Int(input.step_start), Int(input.step_end); vector_scale=vecscale, vector_ref=vecref, output_dir=output_dir)
        return
    end
    error("Config input must include vel or (prefix, step_start, step_end)")
end

function main()
    if length(ARGS) < 1
        println("Usage:")
        println("  julia visualize_step.jl <vel_file>")
        println("  julia visualize_step.jl <vel_prefix> <step_start> <step_end>")
        println("  julia visualize_step.jl <vel_file> --vecscale 0.05")
        println("  julia visualize_step.jl <vel_prefix> <step_start> <step_end> --vecscale 0.05")
        println("  julia visualize_step.jl <vel_file> --vecref 5.0")
        println("  julia visualize_step.jl --config visualize.json")
        println("")
        println("Example:")
        println("  julia visualize_step.jl verification/backward_step/output/vel_0001000.sph")
        println("  julia visualize_step.jl verification/backward_step/output/vel 1000 2000")
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

    vector_scale, vector_ref, args = parse_vector_scale(args)

    if length(args) == 1
        process_file(args[1]; vector_scale=vector_scale, vector_ref=vector_ref)
    elseif length(args) >= 3
        prefix = args[1]
        s_start = parse(Int, args[2])
        s_end = parse(Int, args[3])
        process_range(prefix, s_start, s_end; vector_scale=vector_scale, vector_ref=vector_ref)
    else
        println("Error: Invalid arguments.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
