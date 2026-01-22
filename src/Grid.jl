module Grid

using ..Common

export GridConfig, GridData, generate_grid

"""
    GridConfig

格子生成のための設定構造体。
"""
struct GridConfig
    Nx::Int               # X方向内部セル数
    Ny::Int               # Y方向内部セル数
    Nz::Int               # Z方向内部セル数
    Lx::Float64           # X方向領域長さ
    Ly::Float64           # Y方向領域長さ
    Lz::Float64           # Z方向領域長さ (Added for uniform Z generation)
    origin::NTuple{3, Float64}  # 領域基点
    z_type::Symbol        # :uniform or :non_uniform
    z_file::String        # non-uniform時のファイルパス
end

"""
    GridData

格子データを保持する構造体。
WENO3の5点ステンシルに対応し、各軸両側に2セルのゴーストセルを持つ。
"""
struct GridData
    # 配列サイズ（WENO3の5点ステンシル対応、各軸両側2セルのゴースト）
    mx::Int               # Nx + 4
    my::Int               # Ny + 4
    mz::Int               # Nz + 4

    # セル幅
    dx::Float64           # Lx / Nx
    dy::Float64           # Ly / Ny
    dz::Vector{Float64}   # Z方向セル幅（mz要素）

    # 座標
    x::Vector{Float64}    # X座標（mx要素）
    y::Vector{Float64}    # Y座標（my要素）
    z_face::Vector{Float64}    # Z界面座標（Nz+5要素、HALO込み）
    z_center::Vector{Float64}  # Zセル中心座標（mz要素）
end

"""
    read_z_grid_file(filepath::String, Nz::Int, origin_z::Float64)

Z方向格子ファイルを読み込む。
"""
function read_z_grid_file(filepath::String, Nz::Int, origin_z::Float64)::Vector{Float64}
    if !isfile(filepath)
        error("Grid file not found: $filepath")
    end
    
    lines = readlines(filepath)
    coords = Dict{Int, Float64}()
    expected = Nz + 5
    count = nothing
    for line in lines
        line = split(line, "#")[1]
        line = strip(line)
        if isempty(line)
            continue
        end
        if isnothing(count)
            try
                count = parse(Int, line)
            catch
                error("Invalid grid count in file: $line")
            end
            continue
        end
        parts = split(line)
        if length(parts) < 2
            error("Invalid grid line: $line")
        end
        idx = try
            parse(Int, parts[1])
        catch
            error("Invalid grid index: $line")
        end
        val = try
            parse(Float64, parts[2])
        catch
            error("Invalid grid coordinate: $line")
        end
        coords[idx] = val
    end

    if isnothing(count)
        error("Grid file is missing point count.")
    end
    if count != expected
        error("Grid points mismatch: expected $(expected), got $(count)")
    end
    if length(coords) != expected
        error("Grid points count mismatch: expected $(expected), got $(length(coords))")
    end
    for idx in 1:expected
        if !haskey(coords, idx)
            error("Missing grid index: $(idx)")
        end
    end
    
    z_vals = [coords[i] for i in 1:expected]
    return z_vals .+ origin_z
end

"""
    generate_grid(config::GridConfig, dim_params::DimensionParams)::GridData

格子を生成し、無次元化を行う。
"""
function generate_grid(config::GridConfig, dim_params::DimensionParams)::GridData
    Nx, Ny, Nz = config.Nx, config.Ny, config.Nz
    Lx, Ly, Lz = config.Lx, config.Ly, config.Lz
    L0 = dim_params.L0

    mx = Nx + 4
    my = Ny + 4
    mz = Nz + 4

    # --- X方向 (等間隔) ---
    dx_dim = Lx / Nx
    x_dim = [(i - 2.5) * dx_dim + config.origin[1] for i in 1:mx]

    # --- Y方向 (等間隔) ---
    dy_dim = Ly / Ny
    y_dim = [(j - 2.5) * dy_dim + config.origin[2] for j in 1:my]

    # --- Z方向 ---
    z_face_dim = Vector{Float64}(undef, Nz + 5)
    
    if config.z_type == :uniform
        dz_dim = Lz / Nz
        # 界面座標: k=1 -> origin - 2*dz (Ghost start)
        # Real domain: k=3 (idx=3) is start of real domain (z=origin) ?
        # Wait, if z_face defines interfaces.
        # Real domain is from z_face[3] to z_face[Nz+3].
        # Nz cells.
        # z_face[3] = origin[3].
        # z_face[Nz+3] = origin[3] + Lz.
        # Ghost cells: 2 before, 2 after.
        # k=1: origin - 2*dz
        # k=2: origin - 1*dz
        # k=3: origin
        # ...
        for k in 1:(Nz+5)
            z_face_dim[k] = config.origin[3] + (k - 3) * dz_dim
        end
    elseif config.z_type == :non_uniform
        z_face_dim = read_z_grid_file(config.z_file, Nz, config.origin[3])
    else
        error("Unknown z_type: $(config.z_type)")
    end

    # --- 無次元化 ---
    dx = dx_dim / L0
    dy = dy_dim / L0
    
    x = x_dim ./ L0
    y = y_dim ./ L0
    z_face = z_face_dim ./ L0
    
    z_center = Vector{Float64}(undef, mz)
    dz = Vector{Float64}(undef, mz)
    
    for k in 1:mz
        dz[k] = z_face[k+1] - z_face[k]
        z_center[k] = (z_face[k] + z_face[k+1]) / 2.0
    end

    return GridData(mx, my, mz, dx, dy, dz, x, y, z_face, z_center)
end

end # module Grid
