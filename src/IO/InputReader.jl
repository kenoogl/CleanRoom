module InputReader

using JSON3
using ..Common
using ..Grid
using ..Fields: estimate_memory_size
using ..TimeIntegration
using ..BoundaryConditions: Outflow, Periodic, Symmetric, Wall, SlidingWall, Inflow, Opening, OpeningFlowType, OpeningInlet, OpeningOutlet, ExternalBC, InternalBoundary, BoundaryConditionSet
using ..PressureSolver

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
    debug::Bool
    start::Symbol         # :initial or :restart
    max_step::Int
    dim_params::DimensionParams
    grid_config::GridConfig
    courant_number::Float64
    intervals::IntervalConfig
    poisson::PoissonConfig
    time_scheme::Symbol   # :Euler, :RK2, :RK4
    smagorinsky_constant::Float64
    div_max_threshold::Float64
    initial_condition::InitialCondition
    restart_file::String
    reverse_flow_stabilization::Bool
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
    if value isa Bool
        return value
    end
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
    debug = parse_yesno(get(data, :debug, get(data, :Debug, "no")), "debug")
    reverse_flow_stabilization = parse_yesno(
        get(data, :reverse_flow_stabilization, get(data, :Reverse_flow_stabilization, "no")),
        "reverse_flow_stabilization"
    )
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

    # --- Weak Compressibility (Mach number) ---
    mach_number = if haskey(data, :Mach_number)
        Float64(data.Mach_number)
    elseif haskey(data, :mach_number)
        Float64(data.mach_number)
    elseif haskey(data, :Sound_Speed)
        c = Float64(data.Sound_Speed)
        if c <= 0.0
            error("Sound_Speed must be > 0.")
        end
        U0 / c
    elseif haskey(data, :sound_speed)
        c = Float64(data.sound_speed)
        if c <= 0.0
            error("sound_speed must be > 0.")
        end
        U0 / c
    else
        0.0
    end
    if mach_number < 0.0
        error("Mach_number must be >= 0.")
    end
    mach2 = mach_number^2

    # --- Grid ---
    origin = parse_tuple3(data.Origin_of_Region)
    dom = data.Domain
    z_grid = data.Z_grid
    z_type_str = lowercase(String(z_grid.type))
    z_type = if z_type_str == "uniform"
        :uniform
    elseif z_type_str == "non-uniform" || z_type_str == "non_uniform"
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
    time_scheme_str = lowercase(String(get(data, :Time_Integration_Scheme, "euler")))
    time_scheme = if time_scheme_str == "euler"; :Euler
                  elseif time_scheme_str == "rk2"; :RK2
                  elseif time_scheme_str == "rk4"; :RK4
                  else; error("Unknown Time_Integration_Scheme: $(time_scheme_str)")
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
    if haskey(data, :Visualization)
        @warn "Visualization section is deprecated in solver JSON; use external visualization config instead."
    end

    # --- Poisson ---
    p = data.Poisson_parameter
    p_solver_str = lowercase(String(get(p, :solver, "redbacksor")))
    solver_type = if p_solver_str == "cg"; CG
                  elseif p_solver_str == "bicgstab"; BiCGSTAB
                  else; RedBlackSOR
                  end
    
    div_action_str = lowercase(String(get(p, :on_divergence, "warncontinue")))
    action_type = (div_action_str == "abort") ? Abort : WarnContinue

    precond_str = lowercase(String(get(p, :preconditioner, "sor")))
    precond_type = if precond_str == "none"; PrecondNone
                   elseif precond_str == "sor"; PrecondSOR
                   else; error("Unknown preconditioner: $(precond_str)")
                   end

    poisson_config = PoissonConfig(
        solver_type,
        Float64(p.coef_acceleration),
        Float64(p.convergence_criteria),
        Int(p.Iteration_max),
        action_type,
        precond_type,
        mach2
    )

    # --- Initial Condition ---
    ic = get(data, :Initial_Condition, nothing)
    initial_condition = if isnothing(ic)
        InitialCondition((0.0, 0.0, 0.0), 0.0)
    else
        if haskey(ic, :pressure)
            @warn "Initial_Condition.pressure is deprecated; default p=0 is recommended (pressure will be read if specified)."
        end
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
        dry_run, debug, start_mode, max_step, dim_params,
        grid_config, courant_number, intervals,
        poisson_config, time_scheme,
        smagorinsky_constant, div_max_threshold,
        initial_condition, restart_file, reverse_flow_stabilization
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
        if bc_data isa AbstractString
            error("external_boundaries must be an object with velocity.type/value; string shorthand is not allowed.")
        end
        if haskey(bc_data, :type)
            error("external_boundaries must use velocity.type (top-level \"type\" is not allowed).")
        end
        if !haskey(bc_data, :velocity)
            error("external_boundaries must specify \"velocity\" object.")
        end

        vel_obj = bc_data.velocity
        if vel_obj isa AbstractString
            error("external_boundaries.velocity must be an object with type/value.")
        end
        if !haskey(vel_obj, :type)
            error("external_boundaries.velocity must specify \"type\".")
        end

        type_str = lowercase(String(vel_obj.type))
        # 廃止キーワードのチェック
        if type_str in ["inlet", "outlet", "dirichlet", "neumann"]
            error("Deprecated BC keyword '$(type_str)'. Use 'inflow'/'outflow' instead. For openings, use 'openings' array.")
        end
        # ExternalBCでのOpening禁止
        if type_str == "opening"
            error("'Opening' cannot be specified in external_boundaries. Use 'openings' array instead.")
        end
        
        bc_type = if type_str == "inflow"; Inflow
                  elseif type_str == "outflow"; Outflow
                  elseif type_str == "periodic"; Periodic
                  elseif type_str == "symmetric"; Symmetric
                  elseif type_str == "wall"; Wall
                  elseif type_str == "slidingwall" || type_str == "sliding_wall"; SlidingWall
                  else; error("Unknown velocity BC type: $(type_str)")
                  end
        val = (0.0, 0.0, 0.0)
        if bc_type == Inflow || bc_type == SlidingWall
            if !haskey(vel_obj, :value)
                error("$(bc_type) boundary requires velocity.value field.")
            end
            val_dim = parse_tuple3(vel_obj.value)
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

    openings = Opening[]
    if haskey(data, :openings)
        for item in data.openings
            name = String(get(item, :name, ""))
            shape = Symbol(lowercase(String(item.type)))
            boundary = Symbol(lowercase(String(item.boundary)))
            flow_str = lowercase(String(item.flow_type))
            flow_type = (flow_str == "inlet") ? OpeningInlet : OpeningOutlet

            pos = haskey(item, :position) ? parse_tuple3(item.position) ./ L0 : (0.0, 0.0, 0.0)
            sz = haskey(item, :size) ? parse_tuple2(item.size) ./ L0 : (0.0, 0.0)
            center = haskey(item, :center) ? parse_tuple3(item.center) ./ L0 : (0.0, 0.0, 0.0)
            radius = haskey(item, :radius) ? Float64(item.radius) / L0 : 0.0

            vel = if haskey(item, :velocity)
                parse_tuple3(item.velocity) ./ U0
            else
                (NaN, NaN, NaN)
            end
            if flow_type == OpeningInlet && any(isnan, vel)
                error("Opening inlet requires velocity.")
            end
            push!(openings, Opening(name, shape, boundary, pos, sz, center, radius, flow_type, vel))
        end
    end

    internal_bcs = InternalBoundary[]
    if haskey(data, :internal_boundaries)
        for item in data.internal_boundaries
            shape_str = lowercase(String(item.type))
            vel = parse_tuple3(item.velocity) ./ U0
            if shape_str == "rectangular"
                region = item.region
                rmin = parse_tuple3(region.min) ./ L0
                rmax = parse_tuple3(region.max) ./ L0
                push!(internal_bcs, InternalBoundary(
                    :rectangular, rmin, rmax,
                    (0.0,0.0,0.0), 0.0, 0.0, :z,
                    vel
                ))
            elseif shape_str == "cylindrical"
                cent = parse_tuple3(item.center) ./ L0
                rad = Float64(item.radius) / L0
                h = Float64(item.height) / L0
                ax_str = lowercase(String(item.axis))
                ax = if ax_str == "x"; :x
                     elseif ax_str == "y"; :y
                     else; :z
                     end
                push!(internal_bcs, InternalBoundary(
                    :cylindrical, (0.0,0.0,0.0), (0.0,0.0,0.0),
                    cent, rad, h, ax,
                    vel
                ))
            else
                error("Unknown internal boundary type: $(shape_str)")
            end
        end
    end

    return BoundaryConditionSet(
        x_min, x_max, y_min, y_max, z_min, z_max,
        openings, internal_bcs
    )
end

end # module InputReader
