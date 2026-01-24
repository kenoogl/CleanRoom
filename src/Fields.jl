module Fields

using ..Common
using ..Grid

export CFDBuffers, KrylovBuffers, RKBuffers
export update_time_average!, estimate_memory_size, check_memory_availability, allocate_rk_buffers

"""
    CFDBuffers

物理量配列、作業配列、統計量を一括管理する構造体。
"""
struct CFDBuffers
    # 速度場
    u::Array{Float64, 3}      # u成分
    v::Array{Float64, 3}      # v成分
    w::Array{Float64, 3}      # w成分

    # 圧力場
    p::Array{Float64, 3}      # 圧力
    p_prev::Array{Float64, 3} # 前ステップ圧力（弱圧縮性項用）

    # 擬似速度（Fractional Step用）
    u_star::Array{Float64, 3}
    v_star::Array{Float64, 3}
    w_star::Array{Float64, 3}

    # 時間平均（Welford法によるインクリメンタル平均）
    u_avg::Array{Float64, 3}
    v_avg::Array{Float64, 3}
    w_avg::Array{Float64, 3}
    avg_count::Base.RefValue{Int}   # 平均化サンプル数（n） (MutableにするためRefValue)

    # マスク・物性
    mask::Array{Float64, 3}   # 流体=1, 物体=0

    # ポアソンソルバー作業配列（共通）
    rhs::Array{Float64, 3}    # ポアソン右辺

    # フラックス作業配列（対流+拡散を共用）
    flux_u::Array{Float64, 3}  # 総フラックスu成分（対流+拡散）
    flux_v::Array{Float64, 3}  # 総フラックスv成分（対流+拡散）
    flux_w::Array{Float64, 3}  # 総フラックスw成分（対流+拡散）

    # 乱流粘性
    nu_t::Array{Float64, 3}    # 乱流粘性係数
    nu_eff::Array{Float64, 3}  # 有効粘性係数（分子+乱流）

    # セルフェイス速度（発散計算用）
    u_face_x::Array{Float64, 3}  # X面上のu速度
    v_face_y::Array{Float64, 3}  # Y面上のv速度
    w_face_z::Array{Float64, 3}  # Z面上のw速度
end

"""
    Base.show(io::IO, b::CFDBuffers)

簡易表示。
"""
function Base.show(io::IO, b::CFDBuffers)
    dims = size(b.u)
    print(io, "CFDBuffers(size=$(dims))")
end

"""
    KrylovBuffers

CG法/BiCGSTAB法用の追加作業配列。
"""
struct KrylovBuffers
    r::Array{Float64, 3}      # 残差
    p::Array{Float64, 3}      # 探索方向
    q::Array{Float64, 3}      # 作業配列（Ap）
    # BiCGSTAB追加
    r0::Array{Float64, 3}     # 初期残差（BiCGSTABのみ）
    s::Array{Float64, 3}      # 作業配列（BiCGSTABのみ）
    t::Array{Float64, 3}      # 作業配列（BiCGSTABのみ）
end

"""
    RKBuffers

Runge-Kutta法用の作業配列。
"""
struct RKBuffers
    # RK2: u_rk1のみ使用（初期値保存用）
    # RK4: u_rk1〜u_rk4を使用（初期値+勾配保存用）
    u_rk1::Array{Float64, 3}
    v_rk1::Array{Float64, 3}
    w_rk1::Array{Float64, 3}
    u_rk2::Array{Float64, 3}  # RK4のみ
    v_rk2::Array{Float64, 3}
    w_rk2::Array{Float64, 3}
    u_rk3::Array{Float64, 3}  # RK4のみ
    v_rk3::Array{Float64, 3}
    w_rk3::Array{Float64, 3}
    u_rk4::Array{Float64, 3}  # RK4のみ
    v_rk4::Array{Float64, 3}
    w_rk4::Array{Float64, 3}
end

"""
    allocate_rk_buffers(scheme, mx, my, mz)

RK法用バッファを確保する。Eulerの場合はnothingを返す。
"""
function allocate_rk_buffers(scheme::Symbol, mx::Int, my::Int, mz::Int)
    zeros_3d() = zeros(Float64, mx, my, mz)
    empty_3d() = Array{Float64}(undef, 0, 0, 0)
    
    if scheme == :RK2
        return RKBuffers(
            zeros_3d(), zeros_3d(), zeros_3d(),
            empty_3d(), empty_3d(), empty_3d(),
            empty_3d(), empty_3d(), empty_3d(),
            empty_3d(), empty_3d(), empty_3d()
        )
    elseif scheme == :RK4
        return RKBuffers(
            zeros_3d(), zeros_3d(), zeros_3d(),
            zeros_3d(), zeros_3d(), zeros_3d(),
            zeros_3d(), zeros_3d(), zeros_3d(),
            zeros_3d(), zeros_3d(), zeros_3d()
        )
    else
        return nothing
    end
end

"""
    CFDBuffers(mx::Int, my::Int, mz::Int)

CFDBuffersコンストラクタ。全配列を0で初期化して確保する。
"""
function CFDBuffers(mx::Int, my::Int, mz::Int)
    zeros_3d() = zeros(Float64, mx, my, mz)
    
    return CFDBuffers(
        zeros_3d(), zeros_3d(), zeros_3d(), # u, v, w
        zeros_3d(),                         # p
        zeros_3d(),                         # p_prev
        zeros_3d(), zeros_3d(), zeros_3d(), # u_star...
        zeros_3d(), zeros_3d(), zeros_3d(), # u_avg...
        Ref(0),                             # avg_count
        ones(Float64, mx, my, mz),          # mask (default to fluid=1)
        zeros_3d(),                         # rhs
        zeros_3d(), zeros_3d(), zeros_3d(), # flux...
        zeros_3d(), zeros_3d(),             # nu_t, nu_eff
        zeros_3d(), zeros_3d(), zeros_3d()  # u_face...
    )
end

"""
    update_time_average!(buffers::CFDBuffers, par::String)

Welford法による時間平均の更新。
"""
function update_time_average!(buffers::CFDBuffers, par::String)
    n = buffers.avg_count[] + 1
    buffers.avg_count[] = n
    
    # helper for clean broadcasting
    # mean_new = mean_old + (x - mean_old) / n
    @inbounds begin
        # Loop fusion handled by Julia compiler ideally, or use explicit loops for FLoops if heavy
        # For now using broadcast as it is simple and fast enough for prototype
        @. buffers.u_avg += (buffers.u - buffers.u_avg) / n
        @. buffers.v_avg += (buffers.v - buffers.v_avg) / n
        @. buffers.w_avg += (buffers.w - buffers.w_avg) / n
    end
end

"""
    estimate_memory_size(Nx, Ny, Nz; time_scheme=:Euler, solver=:RedBlackSOR)

メモリ使用量を見積もる（Byte単位）。
"""
function estimate_memory_size(
    Nx::Int, Ny::Int, Nz::Int;
    time_scheme::Symbol = :Euler,
    solver::Symbol = :RedBlackSOR
)::NamedTuple{(:bytes, :gb, :arrays, :breakdown), Tuple{Int, Float64, Int, String}}
    mx, my, mz = Nx + 4, Ny + 4, Nz + 4
    cell_count = mx * my * mz
    bytes_per_cell = 8  # Float64

    # CFDBuffers（必須）: 21配列（p_prev含む）
    base_arrays = 21

    # KrylovBuffers（反復法依存）
    krylov_arrays = 0
    if solver == :CG
        krylov_arrays = 3
    elseif solver == :BiCGSTAB
        krylov_arrays = 6
    end

    # RKBuffers（時間積分スキーム依存）
    rk_arrays = 0
    if time_scheme == :RK2
        rk_arrays = 3
    elseif time_scheme == :RK4
        rk_arrays = 12
    end

    num_arrays = base_arrays + krylov_arrays + rk_arrays
    total_bytes = cell_count * bytes_per_cell * num_arrays
    total_gb = total_bytes / (1024^3)

    breakdown = "CFDBuffers: $(base_arrays), Krylov: $(krylov_arrays), RK: $(rk_arrays)"

    return (bytes = total_bytes, gb = total_gb, arrays = num_arrays, breakdown = breakdown)
end

"""
    check_memory_availability(Nx::Int, Ny::Int, Nz::Int)

メモリが足りるかチェックする（簡易実装）。
"""
function check_memory_availability(Nx::Int, Ny::Int, Nz::Int)::Tuple{Bool, String}
    # 実装: Sys.free_memory() と比較
    # ここでは仮実装として常にTrueを返すが、本来はSys.free_memory()を使う
    estimate = estimate_memory_size(Nx, Ny, Nz)
    req_gb = estimate.gb
    
    # Sys.free_memory() is available in Julia
    free_bytes = Sys.free_memory()
    free_gb = free_bytes / (1024^3)
    
    # 80% usage warning logic could be here
    if req_gb > free_gb * 0.8
        return (true, "Warning: $(round(req_gb, digits=2)) GB required, which is > 80% of free memory ($(round(free_gb, digits=2)) GB)")
    elseif req_gb > free_gb
         return (false, "Insufficient memory: $(round(req_gb, digits=2)) GB required, $(round(free_gb, digits=2)) GB free")
    end
    
    return (true, "$(round(req_gb, digits=2)) GB required")
end

end # module Fields
