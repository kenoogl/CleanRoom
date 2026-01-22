module InputReader

using JSON3
using ..Common
using ..Grid
using ..Fields: estimate_memory_size
using ..TimeIntegration
using ..BoundaryConditions
using ..PressureSolver
using ..Visualization

export IntervalConfig, InitialCondition, SimulationParams
export load_parameters, load_boundary_conditions

"""
    IntervalConfig

出力間隔設定。
"""
struct IntervalConfig
    display::Int              # コンソール表示間隔
    history::Int              # history出力間隔
    instantaneous::Int        # 瞬時値ファイル出力間隔
    averaged::Int             # 平均値ファイル出力間隔
    checkpoint::Int           # チェックポイント間隔
    average_start_time_nd::Float64  # 平均化開始時刻（無次元）
end

"""
    InitialCondition

初期条件設定。
"""
struct InitialCondition
    velocity::NTuple{3, Float64}  # 初期速度ベクトル [m/s]
    pressure::Float64             # 初期圧力 [Pa]
end

"""
    SimulationParams

全シミュレーションパラメータ。
"""
struct SimulationParams
    dry_run::Bool
    start::Symbol         # :initial or :restart
    max_step::Int
    dim_params::DimensionParams
    grid_config::GridConfig
    courant_number::Float64
    intervals::IntervalConfig
    visualization::VizConfig
    poisson::PoissonConfig
    time_scheme::Symbol   # :Euler, :RK2, :RK4
    output_dimensional::Bool
    initial_condition::InitialCondition
    restart_file::String
end

# --- Helper functions for JSON parsing ---

function parse_tuple3(arr)
    return (Float64(arr[1]), Float64(arr[2]), Float64(arr[3]))
end

function parse_tuple3_int(arr)
    return (Int(arr[1]), Int(arr[2]), Int(arr[3]))
end

function parse_tuple2(arr)
    return (Float64(arr[1]), Float64(arr[2]))
end

"""
    load_parameters(filepath::String)::SimulationParams

JSONファイルからパラメータを読み込む。
"""
function load_parameters(filepath::String)::SimulationParams
    if !isfile(filepath)
        error("Parameter file not found: $filepath")
    end
    
    json_str = read(filepath, String)
    data = JSON3.read(json_str)

    # --- Basic Params ---
    dry_run = get(data, :dry_run, false)
    start_mode = Symbol(get(data, :start, "initial"))
    max_step = data.max_step
    restart_file = get(data, :restart_file, "")
    output_dimensional = get(data, :output_dimensional, true)

    # --- Dimensions ---
    dims = data.Dimensions
    L0 = Float64(dims.L0)
    U0 = Float64(dims.U0)
    nu = Float64(dims.nu)
    Re = U0 * L0 / nu
    T0 = L0 / U0
    dim_params = DimensionParams(L0, U0, nu, Re, T0)

    # --- Grid ---
    g = data.Grid
    origin = parse_tuple3(g.origin)
    z_type = Symbol(g.z_type)
    z_file = get(g, :z_file, "")
    # Spec update: Lz is needed if uniform
    Lz = haskey(g, :Lz) ? Float64(g.Lz) : 0.0 # Handle optional if non-uniform?
    
    grid_config = GridConfig(
        Int(g.Nx), Int(g.Ny), Int(g.Nz),
        Float64(g.Lx), Float64(g.Ly), Lz,
        origin, z_type, z_file
    )

    # --- Time & Stability ---
    courant_number = Float64(data.Time.Co)
    time_scheme = Symbol(get(data.Time, :scheme, "Euler"))

    # --- Intervals ---
    intv = data.Intervals
    avg_start_time = Float64(intv.Start_time_for_averaging)
    avg_start_time_nd = avg_start_time / T0
    
    intervals = IntervalConfig(
        Int(intv.display),
        Int(intv.history),
        Int(get(intv, :Instantaneous_file, 0)),
        Int(get(intv, :averaged_file, 0)),
        Int(get(intv, :checkpoint, 0)),
        avg_start_time_nd
    )

    # --- Visualization ---
    viz = get(data, :Visualization, nothing)
    viz_config = if isnothing(viz)
        VizConfig(0, :none, 0, [], :none, "", false, 1, false)
    else
        VizConfig(
            Int(get(viz, :interval, 0)),
            Symbol(get(viz, :plane, "xy")),
            Int(get(viz, :plane_index, 1)),
            Symbol.(get(viz, :variables, ["velocity", "pressure"])),
            Symbol(get(viz, :output_format, "png")),
            String(get(viz, :output_dir, "viz")),
            Bool(get(viz, :vector_enabled, false)),
            Int(get(viz, :vector_skip, 1)),
            Bool(get(viz, :text_output, false))
        )
    end

    # --- Poisson ---
    p = data.Poisson
    p_solver = Symbol(get(p, :solver, "RedBlackSOR"))
    solver_type = if p_solver == :CG; CG
                  elseif p_solver == :BiCGSTAB; BiCGSTAB
                  else; RedBlackSOR
                  end
    
    div_action = Symbol(get(p, :on_divergence, "WarnContinue"))
    action_type = (div_action == :Abort) ? Abort : WarnContinue

    poisson_config = PoissonConfig(
        solver_type,
        Float64(p.omega),
        Float64(p.tol),
        Int(p.max_iter),
        action_type
    )

    # --- Initial Condition ---
    ic = get(data, :InitialCondition, nothing)
    initial_condition = if isnothing(ic)
        InitialCondition((0.0, 0.0, 0.0), 0.0)
    else
        InitialCondition(
            parse_tuple3(get(ic, :velocity, [0.0, 0.0, 0.0])),
            Float64(get(ic, :pressure, 0.0))
        )
    end

    # 安定条件の簡易チェック (D < 0.5)
    # ここでは格子幅が確定していないため完全なチェックはGrid生成後に行うが
    # パラメータ読み込み時点での不整合があればエラーにする

    return SimulationParams(
        dry_run, start_mode, max_step, dim_params,
        grid_config, courant_number, intervals,
        viz_config, poisson_config, time_scheme,
        output_dimensional, initial_condition, restart_file
    )
end

"""
    load_boundary_conditions(filepath::String)::BoundaryConditionSet

境界条件JSONを読み込む。
"""
function load_boundary_conditions(filepath::String)::BoundaryConditionSet
    if !isfile(filepath)
        error("Boundary condition file not found: $filepath")
    end

    json_str = read(filepath, String)
    data = JSON3.read(json_str)

    function parse_bc(bc_data)
        t = Symbol(bc_data.type)
        val = haskey(bc_data, :value) ? parse_tuple3(bc_data.value) : (0.0, 0.0, 0.0)
        
        bc_type = if t == :Dirichlet; Dirichlet
                  elseif t == :Neumann; Neumann
                  elseif t == :Outflow; Outflow
                  elseif t == :Periodic; Periodic
                  else; error("Unknown BC type: $t")
                  end
        return ExternalBC(bc_type, val)
    end

    x_min = parse_bc(data.faces.x_min)
    x_max = parse_bc(data.faces.x_max)
    y_min = parse_bc(data.faces.y_min)
    y_max = parse_bc(data.faces.y_max)
    z_min = parse_bc(data.faces.z_min)
    z_max = parse_bc(data.faces.z_max)

    inlets = InletOutlet[]
    if haskey(data, :Inlets)
        for item in data.Inlets
            push!(inlets, InletOutlet(
                Symbol(item.shape),
                parse_tuple3(item.position),
                parse_tuple2(item.size),
                parse_tuple3_int(item.normal),
                Dirichlet, # TODO: config might specify type? Usually inlets are Dirichlet
                parse_tuple3(item.velocity)
            ))
        end
    end

    outlets = InletOutlet[]
    if haskey(data, :Outlets)
        for item in data.Outlets
            cond = Symbol(get(item, :condition, "Outflow"))
            bctype = (cond == :Dirichlet) ? Dirichlet : Outflow
             # TODO: size can be array of 2 floats or radius?
             # Design doc says "size::NTuple{2, Float64}"
             # If cylinder, maybe size[1]=radius?
             # I'll Assume JSON matches NTuple structure or handle it.
             # If shape is cylindrical, item.radius might exist.
             sz = if Symbol(item.shape) == :cylindrical
                 (Float64(get(item, :radius, 0.0)), 0.0)
             else
                 parse_tuple2(item.size)
             end

            push!(outlets, InletOutlet(
                Symbol(item.shape),
                parse_tuple3(item.position),
                sz,
                parse_tuple3_int(item.normal),
                bctype,
                (0.0,0.0,0.0) # velocity ignored for outflow usually, or 0
            ))
        end
    end

    internal_bcs = InternalBoundary[]
    if haskey(data, :InternalBoundaries)
        for item in data.InternalBoundaries
            shape = Symbol(item.type)
            norm = parse_tuple3_int(item.normal)
            vel = parse_tuple3(item.velocity)
            
            if shape == :rectangular
                rmin = parse_tuple3(item.region_min)
                rmax = parse_tuple3(item.region_max)
                push!(internal_bcs, InternalBoundary(
                    shape, rmin, rmax,
                    (0.0,0.0,0.0), 0.0, 0.0, :z, # dummy for cylinder
                    norm, vel
                ))
            elseif shape == :cylindrical
                cent = parse_tuple3(item.center)
                rad = Float64(item.radius)
                h = Float64(item.height)
                ax = Symbol(item.axis)
                push!(internal_bcs, InternalBoundary(
                    shape, (0.0,0.0,0.0), (0.0,0.0,0.0), # dummy for rect
                    cent, rad, h, ax,
                    norm, vel
                ))
            end
        end
    end

    return BoundaryConditionSet(
        x_min, x_max, y_min, y_max, z_min, z_max,
        inlets, outlets, internal_bcs
    )
end

end # module InputReader
