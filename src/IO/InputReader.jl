module InputReader

using JSON3
using ..Common
using ..Grid
using ..Fields: estimate_memory_size
using ..TimeIntegration
using ..BoundaryConditions: Dirichlet, Neumann, Outflow, Periodic, Symmetric, ExternalBC, InletOutlet, InternalBoundary, BoundaryConditionSet
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
    smagorinsky_constant::Float64
    div_max_threshold::Float64
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

function parse_yesno(value, key::String)::Bool
    if value isa AbstractString
        v = lowercase(String(value))
        if v == "yes"
            return true
        elseif v == "no"
            return false
        end
    end
    error("Invalid value for $(key): expected \"yes\" or \"no\".")
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
    dry_run = parse_yesno(get(data, :dry_run, "no"), "dry_run")
    start_str = lowercase(String(get(data, :start, "initial")))
    start_mode = if start_str == "initial"
        :initial
    elseif start_str == "restart"
        :restart
    else
        error("Invalid start mode: $(start_str)")
    end
    max_step = Int(get(data, :Max_step, 0))
    if max_step <= 0
        error("Max_step must be > 0.")
    end

    # --- Dimensions ---
    L0 = Float64(data.Reference_Length)
    U0 = Float64(data.Reference_Velocity)
    nu = Float64(data.Kinematic_Viscosity)
    Re = U0 * L0 / nu
    T0 = L0 / U0
    dim_params = DimensionParams(L0, U0, nu, Re, T0)

    # --- Model Params ---
    smagorinsky_constant = Float64(get(data, :Smagorinsky_Constant, 0.2))
    div_max_threshold = Float64(get(data, :divMax_threshold, 1.0e-3))

    # --- Grid ---
    origin = parse_tuple3(data.Origin_of_Region)
    dom = data.Domain
    z_grid = data.Z_grid
    z_type_str = String(z_grid.type)
    z_type = if z_type_str == "uniform"
        :uniform
    elseif z_type_str == "non-uniform"
        :non_uniform
    else
        error("Unknown Z_grid.type: $(z_type_str)")
    end
    z_file = haskey(z_grid, :file) ? String(z_grid.file) : ""
    Lz = (z_type == :uniform) ? Float64(z_grid.Lz) : 0.0
    if z_type == :uniform && Lz <= 0.0
        error("Z_grid.Lz must be > 0 for uniform grid.")
    elseif z_type == :non_uniform && isempty(z_file)
        error("Z_grid.file is required for non-uniform grid.")
    end
    grid_config = GridConfig(
        Int(dom.Nx), Int(dom.Ny), Int(dom.Nz),
        Float64(dom.Lx), Float64(dom.Ly), Lz,
        origin, z_type, z_file
    )

    # --- Time & Stability ---
    courant_number = Float64(data.Courant_number)
    time_scheme_str = String(get(data, :Time_Integration_Scheme, "Euler"))
    time_scheme = Symbol(time_scheme_str)
    if !(time_scheme in (:Euler, :RK2, :RK4))
        error("Unknown Time_Integration_Scheme: $(time_scheme_str)")
    end

    # --- Intervals ---
    intv = data.Intervals
    avg_start_time = Float64(get(data, :Start_time_for_averaging, 0.0))
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
    p = data.Poisson_parameter
    p_solver = Symbol(get(p, :solver, "RedBlackSOR"))
    solver_type = if p_solver == :CG; CG
                  elseif p_solver == :BiCGSTAB; BiCGSTAB
                  else; RedBlackSOR
                  end
    
    div_action = Symbol(get(p, :on_divergence, "WarnContinue"))
    action_type = (div_action == :Abort) ? Abort : WarnContinue

    poisson_config = PoissonConfig(
        solver_type,
        Float64(p.coef_acceleration),
        Float64(p.convergence_criteria),
        Int(p.Iteration_max),
        action_type
    )

    # --- Initial Condition ---
    ic = get(data, :Initial_Condition, nothing)
    initial_condition = if isnothing(ic)
        InitialCondition((0.0, 0.0, 0.0), 0.0)
    else
        InitialCondition(
            parse_tuple3(get(ic, :velocity, [0.0, 0.0, 0.0])),
            Float64(get(ic, :pressure, 0.0))
        )
    end

    restart_file = ""
    if haskey(data, :Restart)
        restart_file = String(get(data.Restart, :file, ""))
    end
    if start_mode == :restart && isempty(restart_file)
        error("Restart.file must be set when start=\"restart\".")
    end

    # 安定条件の簡易チェック (D < 0.5)
    # ここでは格子幅が確定していないため完全なチェックはGrid生成後に行うが
    # パラメータ読み込み時点での不整合があればエラーにする

    return SimulationParams(
        dry_run, start_mode, max_step, dim_params,
        grid_config, courant_number, intervals,
        viz_config, poisson_config, time_scheme,
        smagorinsky_constant, div_max_threshold,
        initial_condition, restart_file
    )
end

"""
    load_boundary_conditions(filepath::String, dim_params::DimensionParams)::BoundaryConditionSet

境界条件JSONを読み込む。
"""
function load_boundary_conditions(filepath::String, dim_params::DimensionParams)::BoundaryConditionSet
    if !isfile(filepath)
        error("Boundary condition file not found: $filepath")
    end

    L0 = dim_params.L0
    U0 = dim_params.U0

    json_str = read(filepath, String)
    data = JSON3.read(json_str)

    function parse_bc(bc_data)
        vel_str = lowercase(String(bc_data.velocity))
        bc_type = if vel_str == "dirichlet"; Dirichlet
                  elseif vel_str == "neumann"; Neumann
                  elseif vel_str == "outflow"; Outflow
                  elseif vel_str == "periodic"; Periodic
                  elseif vel_str == "symmetric"; Symmetric
                  else; error("Unknown velocity BC type: $(vel_str)")
                  end
        val = (0.0, 0.0, 0.0)
        if bc_type == Dirichlet
            if !haskey(bc_data, :value)
                error("Dirichlet boundary requires value.")
            end
            val_dim = parse_tuple3(bc_data.value)
            val = val_dim ./ U0
        end
        return ExternalBC(bc_type, val)
    end

    x_min = parse_bc(data.external_boundaries.x_min)
    x_max = parse_bc(data.external_boundaries.x_max)
    y_min = parse_bc(data.external_boundaries.y_min)
    y_max = parse_bc(data.external_boundaries.y_max)
    z_min = parse_bc(data.external_boundaries.z_min)
    z_max = parse_bc(data.external_boundaries.z_max)

    inlets = InletOutlet[]
    if haskey(data, :inlets)
        for item in data.inlets
            shape = Symbol(item.type)
            if shape != :rectangular
                error("Unsupported inlet type: $(shape)")
            end
            pos = parse_tuple3(item.position) ./ L0
            sz = parse_tuple2(item.size) ./ L0
            vel = parse_tuple3(item.velocity) ./ U0
            push!(inlets, InletOutlet(
                shape,
                pos,
                sz,
                parse_tuple3_int(item.normal),
                Dirichlet,
                vel
            ))
        end
    end

    outlets = InletOutlet[]
    if haskey(data, :outlets)
        for item in data.outlets
            shape = Symbol(item.type)
            if shape != :rectangular
                error("Unsupported outlet type: $(shape)")
            end
            cond_str = lowercase(String(get(item, :condition, "outflow")))
            bctype = (cond_str == "dirichlet") ? Dirichlet : Outflow
            pos = parse_tuple3(item.position) ./ L0
            sz = parse_tuple2(item.size) ./ L0
            vel = (0.0, 0.0, 0.0)
            if bctype == Dirichlet
                if !haskey(item, :velocity)
                    error("Outlet dirichlet requires velocity.")
                end
                vel = parse_tuple3(item.velocity) ./ U0
            end
            push!(outlets, InletOutlet(
                shape,
                pos,
                sz,
                parse_tuple3_int(item.normal),
                bctype,
                vel
            ))
        end
    end

    internal_bcs = InternalBoundary[]
    if haskey(data, :internal_boundaries)
        for item in data.internal_boundaries
            shape = Symbol(item.type)
            norm = parse_tuple3_int(item.normal)
            vel = parse_tuple3(item.velocity) ./ U0
            if shape == :rectangular
                region = item.region
                rmin = parse_tuple3(region.min) ./ L0
                rmax = parse_tuple3(region.max) ./ L0
                push!(internal_bcs, InternalBoundary(
                    shape, rmin, rmax,
                    (0.0,0.0,0.0), 0.0, 0.0, :z,
                    norm, vel
                ))
            elseif shape == :cylindrical
                cent = parse_tuple3(item.center) ./ L0
                rad = Float64(item.radius) / L0
                h = Float64(item.height) / L0
                ax = Symbol(item.axis)
                push!(internal_bcs, InternalBoundary(
                    shape, (0.0,0.0,0.0), (0.0,0.0,0.0),
                    cent, rad, h, ax,
                    norm, vel
                ))
            else
                error("Unknown internal boundary type: $(shape)")
            end
        end
    end

    return BoundaryConditionSet(
        x_min, x_max, y_min, y_max, z_min, z_max,
        inlets, outlets, internal_bcs
    )
end

end # module InputReader
