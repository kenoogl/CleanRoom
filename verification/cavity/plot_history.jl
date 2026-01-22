using CairoMakie
using DelimitedFiles

function plot_history(history_file, output_file)
    lines_content = readlines(history_file)
    steps = Float64[]
    du = Float64[]
    pres_res = Float64[]
    
    for line in lines_content
        parts = split(line)
        if length(parts) < 7
            continue
        end
        
        # Attempt to parse all parts as Float64
        try
            vals = parse.(Float64, parts)
            # specific check to avoid accidental parsing of numeric-like strings if any
            push!(steps, vals[1])
            push!(du, vals[5])
            push!(pres_res, vals[7])
        catch
            # If parsing fails (e.g. headers), skip the line
            continue
        end
    end
    
    f = Figure(size = (800, 600))
    ax1 = Axis(f[1, 1], title = "Convergence History", ylabel = "dU", yscale = log10)
    lines!(ax1, steps, du, label = "dU", color = :blue)
    axislegend(ax1)
    
    ax2 = Axis(f[2, 1], ylabel = "Pressure Residual", xlabel = "Step", yscale = log10)
    lines!(ax2, steps, pres_res, label = "Pres Res", color = :red)
    axislegend(ax2)
    
    save(output_file, f)
    println("Saved plot to $output_file")
end

plot_history("verification/cavity/output/history.txt", "verification/cavity/output/convergence.png")
