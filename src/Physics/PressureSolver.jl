module PressureSolver

using ..Common
using ..Fields
using ..Grid
using ..BoundaryConditions
using LinearAlgebra # for norm
using FLoops

export SolverType, DivergenceAction, PreconditionerType, PoissonConfig
export RedBlackSOR, CG, BiCGSTAB, WarnContinue, Abort
export PrecondNone, PrecondSOR
export solve_poisson!

const PRECONDITIONER_SWEEPS = 4

@inline function apply_periodic_pressure_if_needed!(
    p::Array{Float64, 3},
    grid::GridData,
    bc_set
)
    if isnothing(bc_set)
        return
    end
    if bc_set.x_min.velocity_type == Periodic && bc_set.x_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :x)
    end
    if bc_set.y_min.velocity_type == Periodic && bc_set.y_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :y)
    end
    if bc_set.z_min.velocity_type == Periodic && bc_set.z_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :z)
    end
end

@inline function remove_rhs_mean_if_singular!(
    rhs::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    alpha::Float64
)::Float64
    if alpha != 0.0
        return 0.0
    end
    sum_b = 0.0
    count = 0.0
    @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2, i in 3:grid.mx-2
        m0 = mask[i, j, k]
        sum_b += rhs[i, j, k] * m0
        count += m0
    end
    if count == 0.0
        return 0.0
    end
    avg_b = sum_b / count
    @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2, i in 3:grid.mx-2
        rhs[i, j, k] -= avg_b * mask[i, j, k]
    end
    return avg_b
end

@inline function restore_rhs_mean!(
    rhs::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    avg_b::Float64
)
    if avg_b == 0.0
        return
    end
    @inbounds for k in 3:grid.mz-2, j in 3:grid.my-2, i in 3:grid.mx-2
        rhs[i, j, k] += avg_b * mask[i, j, k]
    end
end

@enum SolverType begin
    RedBlackSOR
    CG
    BiCGSTAB
end

@enum DivergenceAction begin
    WarnContinue      # 警告を出力して継続
    Abort             # 計算を停止
end

@enum PreconditionerType begin
    PrecondNone
    PrecondSOR
end

struct PoissonConfig
    solver::SolverType
    omega::Float64             # SOR加速係数
    tol::Float64               # 収束判定値
    max_iter::Int              # 最大反復回数
    on_divergence::DivergenceAction  # 収束失敗時の動作
    preconditioner::PreconditionerType  # 前処理種別（CG/BiCGSTABのみ）
    mach2::Float64             # 弱圧縮性係数 M^2
end



"""
    solve_poisson!(buffers, grid, config, bc_set, par)

圧力ポアソン方程式を解く。
"""
function solve_poisson!(
    buffers::CFDBuffers,
    grid::GridData,
    config::PoissonConfig,
    bc_set,  # Union{Nothing, BoundaryConditionSet} - for periodic BC
    par::String,
    alpha::Float64 = 0.0
)::Tuple{Bool, Int, Float64}
    # Solver dispatch
    result = if config.solver == RedBlackSOR
        solve_poisson_sor!(buffers, grid, config, bc_set, par, alpha)
    elseif config.solver == CG
        solve_poisson_cg!(buffers, grid, config, bc_set, par, alpha)
    elseif config.solver == BiCGSTAB
        solve_poisson_bicgstab!(buffers, grid, config, bc_set, par, alpha)
    else
        error("Unknown solver type: $(config.solver)")
    end
    
    converged, iter, residual = result

    # Mean pressure subtraction
    # Since we use Neumann BCs for pressure on all non-periodic boundaries (via mask=0),
    # the pressure solution is unique only up to a constant.
    # We must enforce mean(p) = 0 (or some reference) to prevent drift.
    mx, my, mz = grid.mx, grid.my, grid.mz
    mask = buffers.mask
    p = buffers.p
    
    sum_p = 0.0
    count = 0.0
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        sum_p += p[i, j, k] * mask[i, j, k]
        count += mask[i, j, k]
    end
    
    avg_p = sum_p / count
    @inbounds for k in 1:mz, j in 1:my, i in 1:mx
        m0 = mask[i, j, k]
        p[i, j, k] = m0 * (p[i, j, k] - avg_p) + (1.0-m0)*avg_p
    end

    if !converged && config.on_divergence == WarnContinue
        # Warning printed by caller or Monitor?
    elseif !converged && config.on_divergence == Abort
        error("Poisson solver failed to converge with solver $(config.solver)")
    end

    return result
end

"""
    solve_poisson_sor!(buffers, grid, config, bc_set, par, alpha)
    
Red-Black SOR法による実装
参考: /Users/Daily/Development/H2/src/NonUniform.jl (rbsor_core!, rbsor!, solveSOR!)
"""
function solve_poisson_sor!(
    buffers::CFDBuffers,
    grid::GridData,
    config::PoissonConfig,
    bc_set,
    par::String,
    alpha::Float64
)::Tuple{Bool, Int, Float64}
    
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    p = buffers.p
    rhs = buffers.rhs
    mask = buffers.mask

    iter = 0
    residual = 0.0
    converged = false

    # Ensure BCs are consistent before computing initial residual
    if !isnothing(bc_set)
        apply_pressure_bcs!(p, grid, mask, bc_set)
    end

    res0 = compute_residual_sor(p, rhs, mask, grid, config.omega, alpha)
    if res0 == 0.0
        res0 = 1.0
    end
    
    dx = grid.dx
    dy = grid.dy
    # 非等方係数のため、それぞれの面積・距離成分を計算
    
    # ワーク配列なしで、ループ内で係数を計算する方式（NonUniform.jlのrbsor_core!と同様）
    # これによりメモリ使用量を抑えつつ、キャッシュ効率を上げる
    
    while iter < config.max_iter
        iter += 1
        
        # SOR Update Loop (Red/Black)
        for color in 0:1
            @floop for k in 3:mz-2, j in 3:my-2
                # k-dependent geometry properties
                dz_k_val = grid.dz[k]
                dZ_p = grid.z_center[k+1] - grid.z_center[k]
                dZ_m = grid.z_center[k] - grid.z_center[k-1]
                
                # Z coefficients
                base_cz_p = (dx * dy) / dZ_p
                base_cz_m = (dx * dy) / dZ_m
                base_cy = (dx * dz_k_val) / dy
                base_cx = (dy * dz_k_val) / dx
                vol = dx * dy * dz_k_val

                @simd for i in (3 + (j + k + color + 1) % 2):2:mx-2
                    m0 = mask[i, j, k]
                    
                    # Neighbor masks
                    m_xm = mask[i-1, j, k]
                    m_xp = mask[i+1, j, k]
                    m_ym = mask[i, j-1, k]
                    m_yp = mask[i, j+1, k]
                    m_zm = mask[i, j, k-1]
                    m_zp = mask[i, j, k+1]
                    
                    # Coefficients
                    cond_xm = base_cx * (m_xm * m0)
                    cond_xp = base_cx * (m_xp * m0)
                    cond_ym = base_cy * (m_ym * m0)
                    cond_yp = base_cy * (m_yp * m0)
                    cond_zm = base_cz_m * (m_zm * m0)
                    cond_zp = base_cz_p * (m_zp * m0)
                    
                    # Diagonal term
                    dd = (1.0 - m0) + (cond_xm + cond_xp + cond_ym + cond_yp + cond_zm + cond_zp) + alpha * vol * m0
                    
                    # RHS term
                    b_val = rhs[i, j, k] * vol
                    
                    # Neighbor sum
                    ss = cond_xm * p[i-1, j, k] + cond_xp * p[i+1, j, k] +
                         cond_ym * p[i, j-1, k] + cond_yp * p[i, j+1, k] +
                         cond_zm * p[i, j, k-1] + cond_zp * p[i, j, k+1]
                    
                    # Update
                    pp = p[i, j, k]
                    dp = ((ss - b_val) / dd - pp) * m0
                    p[i, j, k] = pp + config.omega * dp
                end
            end
        end
        
        # Apply BCs immediately after update to ensure ghost cells are consistent
        if !isnothing(bc_set)
            apply_pressure_bcs!(p, grid, mask, bc_set)
        end
        
        # Compute residual with consistent ghost cells
        # Note: omega is not used in residual calculation (just normalization scale if needed, but we use true residual)
        # We reuse compute_residual_sor which calculates true residual norm
        current_res_norm = compute_residual_sor(p, rhs, mask, grid, config.omega, alpha)
        residual = current_res_norm / res0
        
        if residual < config.tol
            converged = true
            break
        end
    end
    
    return (converged, iter, residual)
end

"""
    solve_poisson_cg!(buffers, grid, config, bc_set, par, alpha)
    
Conjugate Gradient法による実装
参考: /Users/Daily/Development/H2/src/NonUniform.jl (CG!)
"""
function solve_poisson_cg!(
    buffers::CFDBuffers,
    grid::GridData,
    config::PoissonConfig,
    bc_set,
    par::String,
    helm_alpha::Float64
)::Tuple{Bool, Int, Float64}
    
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    
    p = buffers.p
    rhs = buffers.rhs
    mask = buffers.mask

    rhs_avg = remove_rhs_mean_if_singular!(rhs, mask, grid, helm_alpha)
    
    # Reuse buffers for CG work vectors
    r = buffers.flux_u       # Residual
    pk = buffers.flux_v      # Search direction
    q = buffers.flux_w       # A * p (work vector)
    z = buffers.nu_t         # Preconditioned residual (temporary)
    
    # 1. Compute initial residual: r = b - Ap
    apply_periodic_pressure_if_needed!(p, grid, bc_set)
    res0 = calc_residual_cg!(r, p, rhs, mask, grid, par, helm_alpha)
    if res0 == 0.0
        restore_rhs_mean!(rhs, mask, grid, rhs_avg)
        return (true, 0, 0.0)
    end
    
    # 2. Preconditioner: z = M^-1 r (Gauss-Seidel)
    apply_preconditioner!(z, r, mask, grid, helm_alpha, config.preconditioner, bc_set)
    
    # 3. p = z
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        pk[i, j, k] = z[i, j, k]
    end
    
    # 4. rho_old = r·z
    rho_old = dot_product_cg(r, z, mask, grid, par)
    
    converged = false
    iter = 0
    residual = 1.0
    
    for itr in 1:config.max_iter
        iter = itr
        
        # 4. q = A * p
        apply_periodic_pressure_if_needed!(pk, grid, bc_set)
        calc_laplacian_cg!(q, pk, mask, grid, par, helm_alpha)
        
        # 5. alpha = rho_old / (p·q)
        pq = dot_product_cg(pk, q, mask, grid, par)
        
        if abs(pq) < 1e-30
            # Breakdown
            break
        end
        
        alpha_k = rho_old / pq
        
        # 6 & 7. x = x + alpha * p,  r = r - alpha * q
        @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k]
            p[i, j, k] += alpha_k * pk[i, j, k] * m0
            r[i, j, k] -= alpha_k * q[i, j, k] * m0
        end
        
        # 8. res = ||r|| / res0
        res_sq = dot_self_cg(r, mask, grid, par)
        residual = sqrt(res_sq) / res0
        
        if residual < config.tol
            converged = true
            break
        end
        
        # 9. Preconditioner: z = M^-1 r
        apply_preconditioner!(z, r, mask, grid, helm_alpha, config.preconditioner, bc_set)
        
        # 10. rho_new = r·z
        rho_new = dot_product_cg(r, z, mask, grid, par)
        
        if abs(rho_new) < 1e-30
            break
        end
        
        # 11. beta = rho_new / rho_old
        beta = rho_new / rho_old
        
        # 12. p = z + beta * p
        @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
            pk[i, j, k] = z[i, j, k] + beta * pk[i, j, k]
        end
        
        rho_old = rho_new
    end
    
    if !isnothing(bc_set)
        apply_pressure_bcs!(p, grid, mask, bc_set)
    end

    restore_rhs_mean!(rhs, mask, grid, rhs_avg)
    
    return (converged, iter, residual)
end


"""
    calc_laplacian_cg!(Ap, p, mask, grid, par, alpha)
    
Compute Ap = +∇²p + α p (SPD operator for weak compressibility)
"""
function calc_laplacian_cg!(
    Ap::Array{Float64, 3},
    p::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    par::String,
    alpha::Float64
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    
    @inbounds for k in 3:mz-2
        dz_k = grid.dz[k]
        dZ_p = grid.z_center[k+1] - grid.z_center[k]
        dZ_m = grid.z_center[k] - grid.z_center[k-1]
        
        base_cz_p = (dx * dy) / dZ_p
        base_cz_m = (dx * dy) / dZ_m
        base_cy = (dx * dz_k) / dy
        base_cx = (dy * dz_k) / dx
        
        vol = dx * dy * dz_k
        for j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k]
            
            m_xm = mask[i-1, j, k]
            m_xp = mask[i+1, j, k]
            m_ym = mask[i, j-1, k]
            m_yp = mask[i, j+1, k]
            m_zm = mask[i, j, k-1]
            m_zp = mask[i, j, k+1]
            
            cond_xm = base_cx * (m_xm * m0)
            cond_xp = base_cx * (m_xp * m0)
            cond_ym = base_cy * (m_ym * m0)
            cond_yp = base_cy * (m_yp * m0)
            cond_zm = base_cz_m * (m_zm * m0)
            cond_zp = base_cz_p * (m_zp * m0)
            
            dd = (1.0 - m0) + (cond_xm + cond_xp + cond_ym + cond_yp + cond_zm + cond_zp) + alpha * vol * m0
            
            ss = cond_xm * p[i-1, j, k] + cond_xp * p[i+1, j, k] +
                 cond_ym * p[i, j-1, k] + cond_yp * p[i, j+1, k] +
                 cond_zm * p[i, j, k-1] + cond_zp * p[i, j, k+1]
            
            # Ap = (-A) * p = (dd * p - ss)  (SPD form)
            Ap[i, j, k] = (dd * p[i, j, k] - ss) * m0
        end
    end
end


"""
    calc_residual_cg!(r, p, rhs, mask, grid, par, alpha)
    
Compute residual r = b - Ap and return ||r|| (b is negated for SPD form)
"""
function calc_residual_cg!(
    r::Array{Float64, 3},
    p::Array{Float64, 3},
    rhs::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    par::String,
    alpha::Float64
)::Float64
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    
    res_sq = 0.0
    
    @inbounds for k in 3:mz-2
        dz_k = grid.dz[k]
        dZ_p = grid.z_center[k+1] - grid.z_center[k]
        dZ_m = grid.z_center[k] - grid.z_center[k-1]
        
        base_cz_p = (dx * dy) / dZ_p
        base_cz_m = (dx * dy) / dZ_m
        base_cy = (dx * dz_k) / dy
        base_cx = (dy * dz_k) / dx
        vol = dx * dy * dz_k
        
        for j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k]
            
            m_xm = mask[i-1, j, k]
            m_xp = mask[i+1, j, k]
            m_ym = mask[i, j-1, k]
            m_yp = mask[i, j+1, k]
            m_zm = mask[i, j, k-1]
            m_zp = mask[i, j, k+1]
            
            cond_xm = base_cx * (m_xm * m0)
            cond_xp = base_cx * (m_xp * m0)
            cond_ym = base_cy * (m_ym * m0)
            cond_yp = base_cy * (m_yp * m0)
            cond_zm = base_cz_m * (m_zm * m0)
            cond_zp = base_cz_p * (m_zp * m0)
            
            dd = (1.0 - m0) + (cond_xm + cond_xp + cond_ym + cond_yp + cond_zm + cond_zp) + alpha * vol * m0
            
            ss = cond_xm * p[i-1, j, k] + cond_xp * p[i+1, j, k] +
                 cond_ym * p[i, j-1, k] + cond_yp * p[i, j+1, k] +
                 cond_zm * p[i, j, k-1] + cond_zp * p[i, j, k+1]
            
            b_val = -rhs[i, j, k] * vol
            
            # Use SPD form: A' = -A, b' = -b
            # A' * p = dd * p - ss, r = b' - A' * p
            r_val = (b_val - (dd * p[i, j, k] - ss)) * m0
            r[i, j, k] = r_val
            res_sq += r_val * r_val
        end
    end
    
    return sqrt(res_sq)
end

"""
    compute_residual_sor(p, rhs, mask, grid, omega, alpha)

Compute SOR residual norm (H2-style) for normalization (Helmholtz対応).
"""
function compute_residual_sor(
    p::Array{Float64, 3},
    rhs::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    omega::Float64,
    alpha::Float64
)::Float64
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    res_sq = 0.0

    @floop for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        dz_k_val = grid.dz[k]
        dZ_p = grid.z_center[k+1] - grid.z_center[k]
        dZ_m = grid.z_center[k] - grid.z_center[k-1]

        base_cz_p = (dx * dy) / dZ_p
        base_cz_m = (dx * dy) / dZ_m
        base_cy = (dx * dz_k_val) / dy
        base_cx = (dy * dz_k_val) / dx
        vol = dx * dy * dz_k_val

        m0 = mask[i, j, k]

        cond_xm = base_cx * (mask[i-1, j, k] * m0)
        cond_xp = base_cx * (mask[i+1, j, k] * m0)
        cond_ym = base_cy * (mask[i, j-1, k] * m0)
        cond_yp = base_cy * (mask[i, j+1, k] * m0)
        cond_zm = base_cz_m * (mask[i, j, k-1] * m0)
        cond_zp = base_cz_p * (mask[i, j, k+1] * m0)

        dd = (1.0 - m0) + (cond_xm + cond_xp + cond_ym + cond_yp + cond_zm + cond_zp) + alpha * vol * m0
        ss = cond_xm * p[i-1, j, k] + cond_xp * p[i+1, j, k] +
             cond_ym * p[i, j-1, k] + cond_yp * p[i, j+1, k] +
             cond_zm * p[i, j, k-1] + cond_zp * p[i, j, k+1]

        b_val = rhs[i, j, k] * vol
        # r = (-b_val) - (dd*p - ss) = -b_val - dd*p + ss
        r_val = (-b_val - dd * p[i, j, k] + ss) * m0
        @reduce(res_sq = 0.0 + r_val * r_val)
    end

    return sqrt(res_sq)
end

"""
    apply_preconditioner!(z, r, mask, grid, alpha, precond)

Apply preconditioner to solve M z = r.
"""
function apply_preconditioner!(
    z::Array{Float64, 3},
    r::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    alpha::Float64,
    precond::PreconditionerType,
    bc_set=nothing
)
    if precond == PrecondNone
        copyto!(z, r)
        return
    end
    fill!(z, 0.0)
    dx, dy = grid.dx, grid.dy
    mx, my, mz = grid.mx, grid.my, grid.mz
    omega = 1.0

    for _ in 1:PRECONDITIONER_SWEEPS
        for color in 0:1
            @floop for k in 3:mz-2, j in 3:my-2
                dz_k_val = grid.dz[k]
                dZ_p = grid.z_center[k+1] - grid.z_center[k]
                dZ_m = grid.z_center[k] - grid.z_center[k-1]

                base_cz_p = (dx * dy) / dZ_p
                base_cz_m = (dx * dy) / dZ_m
                base_cy = (dx * dz_k_val) / dy
                base_cx = (dy * dz_k_val) / dx
                vol = dx * dy * dz_k_val
                @simd for i in (3 + (j + k + color + 1) % 2):2:mx-2
                    m0 = mask[i, j, k]
                    
                    cond_xm = base_cx * (mask[i-1, j, k] * m0)
                    cond_xp = base_cx * (mask[i+1, j, k] * m0)
                    cond_ym = base_cy * (mask[i, j-1, k] * m0)
                    cond_yp = base_cy * (mask[i, j+1, k] * m0)
                    cond_zm = base_cz_m * (mask[i, j, k-1] * m0)
                    cond_zp = base_cz_p * (mask[i, j, k+1] * m0)

                    dd = (1.0 - m0) + (cond_xm + cond_xp + cond_ym + cond_yp + cond_zm + cond_zp) + alpha * vol * m0
                    ss = cond_xm * z[i-1, j, k] + cond_xp * z[i+1, j, k] +
                         cond_ym * z[i, j-1, k] + cond_yp * z[i, j+1, k] +
                         cond_zm * z[i, j, k-1] + cond_zp * z[i, j, k+1]

                    b_val = r[i, j, k]
                    dp = ((ss - b_val) / dd - z[i, j, k]) * m0
                    z[i, j, k] += omega * dp
                end
            end
        end
        apply_periodic_pressure_if_needed!(z, grid, bc_set)
    end
end


"""
    dot_self_cg(x, mask, grid, par)
    
Compute ||x||² = Σ x²
"""
function dot_self_cg(
    x::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    par::String
)::Float64
    mx, my, mz = grid.mx, grid.my, grid.mz
    sum_val = 0.0
    
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        if mask[i, j, k] > 0.0
            sum_val += x[i, j, k] * x[i, j, k]
        end
    end
    
    return sum_val
end


"""
    dot_product_cg(x, y, mask, grid, par)
    
Compute x·y = Σ x * y
"""
function dot_product_cg(
    x::Array{Float64, 3},
    y::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    par::String
)::Float64
    mx, my, mz = grid.mx, grid.my, grid.mz
    sum_val = 0.0
    
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        if mask[i, j, k] > 0.0
            sum_val += x[i, j, k] * y[i, j, k]
        end
    end
    
    return sum_val
end


"""
    solve_poisson_bicgstab!(buffers, grid, config, bc_set, par, alpha)
    
BiCGSTAB法による実装（前処理付き）
"""
function solve_poisson_bicgstab!(
    buffers::CFDBuffers,
    grid::GridData,
    config::PoissonConfig,
    bc_set,
    par::String,
    alpha::Float64
)::Tuple{Bool, Int, Float64}
    mx, my, mz = grid.mx, grid.my, grid.mz

    p = buffers.p
    rhs = buffers.rhs
    mask = buffers.mask

    rhs_avg = remove_rhs_mean_if_singular!(rhs, mask, grid, alpha)

    # Work vectors (reuse existing buffers)
    r = buffers.flux_u       # residual / s
    r_hat = buffers.p_prev   # shadow residual
    pvec = buffers.flux_v    # search direction
    v = buffers.flux_w       # A * phat
    phat = buffers.nu_t      # M^-1 * p
    shat = buffers.nu_eff    # M^-1 * s
    t = buffers.rhs          # reuse rhs as temporary for t = A * shat

    # Initial residual r = b - A p
    apply_periodic_pressure_if_needed!(p, grid, bc_set)
    res0 = calc_residual_cg!(r, p, rhs, mask, grid, par, alpha)
    if res0 == 0.0
        restore_rhs_mean!(rhs, mask, grid, rhs_avg)
        return (true, 0, 0.0)
    end

    copyto!(r_hat, r)
    fill!(pvec, 0.0)
    fill!(v, 0.0)

    rho_old = 1.0
    alpha_k = 1.0
    omega_k = 1.0
    converged = false
    iter = 0
    residual = 1.0

    for itr in 1:config.max_iter
        iter = itr

        rho_new = dot_product_cg(r_hat, r, mask, grid, par)
        if abs(rho_new) < 1e-30
            break
        end

        if itr == 1
            @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
                pvec[i, j, k] = r[i, j, k]
            end
        else
            beta = (rho_new / rho_old) * (alpha_k / omega_k)
            @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
                pvec[i, j, k] = r[i, j, k] + beta * (pvec[i, j, k] - omega_k * v[i, j, k])
            end
        end

        # phat = M^-1 p
        apply_preconditioner!(phat, pvec, mask, grid, alpha, config.preconditioner, bc_set)

        # v = A * phat
        apply_periodic_pressure_if_needed!(phat, grid, bc_set)
        calc_laplacian_cg!(v, phat, mask, grid, par, alpha)

        denom = dot_product_cg(r_hat, v, mask, grid, par)
        if abs(denom) < 1e-30
            break
        end
        alpha_k = rho_new / denom

        # s = r - alpha * v (store in r)
        @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k]
            r[i, j, k] = (r[i, j, k] - alpha_k * v[i, j, k]) * m0
        end

        # Check convergence with s
        residual = sqrt(dot_self_cg(r, mask, grid, par)) / res0
        if residual < config.tol
            @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
                m0 = mask[i, j, k]
                p[i, j, k] += alpha_k * phat[i, j, k] * m0
            end
            converged = true
            break
        end

        # shat = M^-1 s
        apply_preconditioner!(shat, r, mask, grid, alpha, config.preconditioner, bc_set)

        # t = A * shat (reuse rhs array)
        apply_periodic_pressure_if_needed!(shat, grid, bc_set)
        calc_laplacian_cg!(t, shat, mask, grid, par, alpha)

        t_dot_s = dot_product_cg(t, r, mask, grid, par)
        t_dot_t = dot_product_cg(t, t, mask, grid, par)
        if abs(t_dot_t) < 1e-30
            break
        end
        omega_k = t_dot_s / t_dot_t

        # x = x + alpha * phat + omega * shat
        @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k]
            p[i, j, k] += (alpha_k * phat[i, j, k] + omega_k * shat[i, j, k]) * m0
        end

        # r = s - omega * t  (s is in r)
        @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k]
            r[i, j, k] = (r[i, j, k] - omega_k * t[i, j, k]) * m0
        end

        residual = sqrt(dot_self_cg(r, mask, grid, par)) / res0
        if residual < config.tol
            converged = true
            break
        end

        if abs(omega_k) < 1e-30
            break
        end

        rho_old = rho_new
    end

    if !isnothing(bc_set)
        apply_pressure_bcs!(p, grid, mask, bc_set)
    end

    restore_rhs_mean!(rhs, mask, grid, rhs_avg)

    return (converged, iter, residual)
end

end # module PressureSolver
