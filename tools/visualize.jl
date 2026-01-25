#!/usr/bin/env julia

module VizDispatcher

using JSON3

const TOOL_DIR = @__DIR__

function usage()
    println("Usage:")
    println("  julia tools/visualize.jl <pattern> [args...]")
    println("  julia tools/visualize.jl --config <visualize.json>")
    println("")
    println("Patterns:")
    println("  cavity         - Cavity flow profiles (velocity SPH only)")
    println("  backward_step  - Backward-facing step (velocity+pressure SPH)")
    println("")
    println("Examples:")
    println("  julia tools/visualize.jl cavity verification/cavity/output/vel_0001000.sph")
    println("  julia tools/visualize.jl cavity verification/cavity/output/vel 1000 2000")
    println("  julia tools/visualize.jl backward_step verification/backward_step/output/vel_0001000.sph")
    println("  julia tools/visualize.jl backward_step verification/backward_step/output/vel 1000 2000")
    println("  julia tools/visualize.jl --config verification/backward_step/visualize.json")
end

function parse_config_arg(args::Vector{String})
    config_path = nothing
    filtered = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--config" || arg == "-c"
            if i == length(args)
                error("Missing value for --config")
            end
            config_path = args[i + 1]
            i += 2
        elseif startswith(arg, "--config=")
            config_path = split(arg, "=", limit=2)[2]
            i += 1
        else
            push!(filtered, arg)
            i += 1
        end
    end
    return config_path, filtered
end

function run_pattern(pattern::String, args::Vector{String})
    function run_tool(path::String)
        # Run in a fresh Julia process so results match direct tool execution.
        base_cmd = Base.julia_cmd()
        exec = copy(base_cmd.exec)
        project = Base.active_project()
        if project !== nothing && !isempty(project)
            push!(exec, "--project=$(dirname(project))")
        end
        push!(exec, path)
        append!(exec, args)
        run(Cmd(exec))
    end

    if pattern == "cavity"
        run_tool(joinpath(TOOL_DIR, "visualize_cavity.jl"))
    elseif pattern == "backward_step"
        run_tool(joinpath(TOOL_DIR, "visualize_step.jl"))
    else
        error("Unknown pattern: $(pattern)")
    end
end

function main()
    args = copy(ARGS)
    if !isempty(args) && endswith(args[1], "visualize.jl")
        args = args[2:end]
    end

    if length(args) == 0 || args[1] in ("-h", "--help")
        usage()
        return
    end

    config_path, args = parse_config_arg(args)
    if config_path !== nothing
        if !isfile(config_path)
            error("Config file not found: $(config_path)")
        end
        data = JSON3.read(read(config_path, String))
        tool = lowercase(String(get(data, :tool, "")))
        if isempty(tool)
            error("Config must include tool")
        end
        run_pattern(tool, ["--config", config_path])
        return
    end

    pattern = lowercase(args[1])
    args = length(args) >= 2 ? args[2:end] : String[]

    if isempty(args)
        println("Error: Missing arguments for pattern $(pattern)")
        usage()
        return
    end

    run_pattern(pattern, args)
end

end # module VizDispatcher

if abspath(PROGRAM_FILE) == @__FILE__
    VizDispatcher.main()
end
