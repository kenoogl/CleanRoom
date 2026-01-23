#!/usr/bin/env julia
using CairoMakie

include("sph_reader.jl")

"""
    plot_step_flow(vel_data, prs_data; output_path=nothing)

Plot velocity vectors and pressure contours for the backward facing step flow.
"""
function plot_step_flow(vd, pd; output_path=nothing)
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
    p = pd.p[:, j, :]
    
    xrange = x[end] - x[1]
    zrange = z[end] - z[1]
    ratio = xrange / zrange
    base_width = 1000
    title_pad = 120
    plot_height = Int(round(2 * base_width / ratio + title_pad))
    fig = Figure(size = (base_width, plot_height), figure_padding = 10)
    supertitle = "Step Flow (Step $(vd.step), t=$(round(vd.time, digits=2)))"
    Label(fig[1, :], supertitle, fontsize = 20)
    
    # 1. Velocity Vectors (Arrows)
    ax1 = Axis(fig[2, 1], title = "Velocity Vectors", xlabel = "x [m]", ylabel = "z [m]",
               aspect = DataAspect(), limits = ((x[1], x[end]), (z[1], z[end])),
               xautolimitmargin = (0.0, 0.0), yautolimitmargin = (0.0, 0.0))
    
    # Downsample for visibility
    step_vec = 4
    arrows2d!(ax1, x[1:step_vec:end], z[1:step_vec:end],
              u[1:step_vec:end, 1:step_vec:end], w[1:step_vec:end, 1:step_vec:end],
              tipwidth = 7.5, tiplength = 7.5, lengthscale = 0.05, color = :black)
    
    # 2. Pressure Contours
    ax2 = Axis(fig[3, 1], title = "Pressure Contours", xlabel = "x [m]", ylabel = "z [m]",
               aspect = DataAspect(), limits = ((x[1], x[end]), (z[1], z[end])),
               xautolimitmargin = (0.0, 0.0), yautolimitmargin = (0.0, 0.0))
    contour!(ax2, x, z, p, levels = 20, color = :black, linewidth = 0.8)

    colsize!(fig.layout, 1, Relative(1.0))
    rowsize!(fig.layout, 1, Auto(0.06))
    rowsize!(fig.layout, 2, Relative(0.47))
    rowsize!(fig.layout, 3, Relative(0.47))
    rowgap!(fig.layout, 10)
    
    # Geometry Overlay
    for ax in (ax1, ax2)
        # Step boundary (approximate based on standard BFS geometry)
        lines!(ax, [-0.5, 0.0, 0.0], [0.3, 0.3, 0.0], color = :blue, linewidth = 2)
        # Solid fill
        poly!(ax, Point2f[(-0.5, 0.0), (0.0, 0.0), (0.0, 0.3), (-0.5, 0.3)],
              color = (:gray, 0.3), strokewidth = 0)
    end
    
    if !isnothing(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end
    
    return fig
end

"""
    process_file(vel_path)

Process a single velocity file (and corresponding pressure file).
"""
function process_file(vel_path::String)
    println("Processing: $vel_path")
    
    # Infer pressure file path
    # Assume naming convention: vel_XXXXXXX.sph -> prs_XXXXXXX.sph
    dir = dirname(vel_path)
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
    
    plot_step_flow(vd, pd, output_path=out_path)
end

"""
    process_range(prefix, step_start, step_end)

Process a range of steps.
"""
function process_range(prefix::String, step_start::Int, step_end::Int)
    for step in step_start:step_end
        # Construct filename. Prefix is likely "verification/backward_step/output/vel"
        # We need to append "_XXXXXXX.sph"
        
        # Check if prefix already includes underscore or not
        # Standard usage: tools/visualize_step.jl verification/backward_step/output/vel 100 200
        
        fname = "$(prefix)_$(lpad(step, 7, '0')).sph"
        if isfile(fname)
            process_file(fname)
        end
    end
end

function main()
    if length(ARGS) < 1
        println("Usage:")
        println("  julia visualize_step.jl <vel_file>")
        println("  julia visualize_step.jl <vel_prefix> <step_start> <step_end>")
        println("")
        println("Example:")
        println("  julia visualize_step.jl verification/backward_step/output/vel_0001000.sph")
        println("  julia visualize_step.jl verification/backward_step/output/vel 1000 2000")
        return
    end
    
    if length(ARGS) == 1
        process_file(ARGS[1])
    elseif length(ARGS) >= 3
        prefix = ARGS[1]
        s_start = parse(Int, ARGS[2])
        s_end = parse(Int, ARGS[3])
        process_range(prefix, s_start, s_end)
    else
        println("Error: Invalid arguments.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
