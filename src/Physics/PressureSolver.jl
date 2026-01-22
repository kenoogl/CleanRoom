module PressureSolver

using ..Common
using ..Fields
using ..Grid
using ..BoundaryConditions
using LinearAlgebra # for norm
using FLoops

export SolverType, DivergenceAction, PoissonConfig
export RedBlackSOR, CG, BiCGSTAB, WarnContinue, Abort
export solve_poisson!

@enum SolverType begin
    RedBlackSOR
    CG
    BiCGSTAB
end

@enum DivergenceAction begin
    WarnContinue      # 警告を出力して継続
    Abort             # 計算を停止
end

struct PoissonConfig
    solver::SolverType
    omega::Float64             # SOR加速係数
    tol::Float64               # 収束判定値
    max_iter::Int              # 最大反復回数
    on_divergence::DivergenceAction  # 収束失敗時の動作
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
    par::String
)::Tuple{Bool, Int, Float64}
    
    # Solver dispatch
    result = if config.solver == RedBlackSOR
        solve_poisson_sor!(buffers, grid, config, bc_set, par)
    elseif config.solver == CG
        solve_poisson_cg!(buffers, grid, config, bc_set, par)
    elseif config.solver == BiCGSTAB
        solve_poisson_bicgstab!(buffers, grid, config, bc_set, par)
    else
        error("Unknown solver type: $(config.solver)")
    end
    
    converged, iter, residual = result

    # Mean pressure subtraction (Pinning) - Common for all solvers
    # Applied after solving to ensure mean pressure is zero
    mx, my, mz = grid.mx, grid.my, grid.mz
    mask = buffers.mask
    p = buffers.p
    
    sum_p = 0.0
    count = 0
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        if mask[i, j, k] > 0
            sum_p += p[i, j, k]
            count += 1
        end
    end
    
    if count > 0
        avg_p = sum_p / count
        @inbounds for k in 1:mz, j in 1:my, i in 1:mx # Subtract everywhere
             p[i, j, k] -= avg_p
        end
    end

    if !converged && config.on_divergence == WarnContinue
        # Warning printed by caller or Monitor?
    elseif !converged && config.on_divergence == Abort
        error("Poisson solver failed to converge with solver $(config.solver)")
    end

    return result
end

"""
    solve_poisson_sor!(buffers, grid, config, bc_set, par)
    
Red-Black SOR法による実装
参考: /Users/Daily/Development/H2/src/NonUniform.jl (rbsor_core!, rbsor!, solveSOR!)
"""
function solve_poisson_sor!(
    buffers::CFDBuffers,
    grid::GridData,
    config::PoissonConfig,
    bc_set,
    par::String
)::Tuple{Bool, Int, Float64}
    
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    p = buffers.p
    rhs = buffers.rhs
    mask = buffers.mask
    
    iter = 0
    residual = 0.0
    converged = false
    
    dx = grid.dx
    dy = grid.dy
    # 非等方係数のため、それぞれの面積・距離成分を計算
    
    # ワーク配列なしで、ループ内で係数を計算する方式（NonUniform.jlのrbsor_core!と同様）
    # これによりメモリ使用量を抑えつつ、キャッシュ効率を上げる
    
    while iter < config.max_iter
        iter += 1
        total_res_sq = 0.0
        
        # SOR Loop (Red/Black)
        for color in 0:1
            @floop for k in 3:mz-2, j in 3:my-2
                # k-dependent geometry properties
                dz_k_val = grid.dz[k]
                dZ_p = grid.z_center[k+1] - grid.z_center[k]
                dZ_m = grid.z_center[k] - grid.z_center[k-1]
                
                # Z coefficients
                base_cz_p = (dx * dy) / dZ_p
                base_cz_m = (dx * dy) / dZ_m
                
                # Y coefficients
                base_cy = (dx * dz_k_val) / dy
                
                # X coefficients
                base_cx = (dy * dz_k_val) / dx
                
                # Volume element
                vol = dx * dy * dz_k_val

                # Red/Black check
                
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
                    dd = (1.0 - m0) + (cond_xm + cond_xp + cond_ym + cond_yp + cond_zm + cond_zp)
                    
                    # RHS term
                    b_val = rhs[i, j, k] * vol
                    
                    # Neighbor sum
                    ss = cond_xm * p[i-1, j, k] + cond_xp * p[i+1, j, k] +
                         cond_ym * p[i, j-1, k] + cond_yp * p[i, j+1, k] +
                         cond_zm * p[i, j, k-1] + cond_zp * p[i, j, k+1]
                    
                    # Update
                    pp = p[i, j, k]
                    dp = ((ss - b_val) / dd - pp) * m0
                    p_new = pp + config.omega * dp
                    
                    # Residual
                    # NonUniform.jl: r = (dd + ω * (cond_xm + cond_ym + cond_zm)) * dp / ω
                    sum_cond_minus = cond_xm + cond_ym + cond_zm
                    r_val = (dd + config.omega * sum_cond_minus) * dp / config.omega
                    
                    res_sq = r_val * r_val
                    @reduce(local_res_sum = 0.0 + res_sq)
                    
                    p[i, j, k] = p_new
                end
            end
            
            total_res_sq += local_res_sum
        end
        
        # Apply periodic BC
        if !isnothing(bc_set)
            apply_periodic_pressure!(p, grid, bc_set)
        end
        
        residual = sqrt(total_res_sq)
        if residual < config.tol
            converged = true
            break
        end
    end
    
    return (converged, iter, residual)
end

"""
    solve_poisson_cg!(buffers, grid, config, bc_set, par)
    
Conjugate Gradient法による実装
参考: /Users/Daily/Development/H2/src/NonUniform.jl (CG!)
"""
function solve_poisson_cg!(
    buffers::CFDBuffers,
    grid::GridData,
    config::PoissonConfig,
    bc_set,
    par::String
)::Tuple{Bool, Int, Float64}
    
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    
    p = buffers.p
    rhs = buffers.rhs
    mask = buffers.mask
    
    # Reuse flux buffers for CG work vectors
    r = buffers.flux_u    # Residual
    pk = buffers.flux_v   # Search direction
    q = buffers.flux_w    # A * p (work vector)
    
    # 1. Compute initial residual: r = b - Ap
    res0 = calc_residual_cg!(r, p, rhs, mask, grid, par)
    
    if res0 < 1e-14
        return (true, 0, res0)
    end
    
    # 2. p = r (no preconditioner)
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        pk[i, j, k] = r[i, j, k]
    end
    
    # 3. rho_old = r·r
    rho_old = dot_self_cg(r, mask, grid, par)
    
    converged = false
    iter = 0
    residual = res0
    
    for itr in 1:config.max_iter
        iter = itr
        
        # 4. q = A * p
        calc_laplacian_cg!(q, pk, mask, grid, par)
        
        # 5. alpha = rho_old / (p·q)
        pq = dot_product_cg(pk, q, mask, grid, par)
        
        if abs(pq) < 1e-30
            # Breakdown
            break
        end
        
        alpha = rho_old / pq
        
        # 6 & 7. x = x + alpha * p,  r = r - alpha * q
        @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k]
            p[i, j, k] += alpha * pk[i, j, k] * m0
            r[i, j, k] -= alpha * q[i, j, k] * m0
        end
        
        # 8. res = ||r|| / res0
        res_sq = dot_self_cg(r, mask, grid, par)
        residual = sqrt(res_sq) / res0
        
        if residual < config.tol
            converged = true
            break
        end
        
        # 9. rho_new = r·r
        rho_new = res_sq
        
        if abs(rho_new) < 1e-30
            break
        end
        
        # 10. beta = rho_new / rho_old
        beta = rho_new / rho_old
        
        # 11. p = r + beta * p
        @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
            pk[i, j, k] = r[i, j, k] + beta * pk[i, j, k]
        end
        
        rho_old = rho_new
    end
    
    # Apply periodic BC
    if !isnothing(bc_set)
        apply_periodic_pressure!(p, grid, bc_set)
    end
    
    return (converged, iter, residual)
end


"""
    calc_laplacian_cg!(Ap, p, mask, grid, par)
    
Compute Ap = -∇²p (Laplacian operator for pressure Poisson equation)
"""
function calc_laplacian_cg!(
    Ap::Array{Float64, 3},
    p::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    par::String
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
            
            dd = (1.0 - m0) + (cond_xm + cond_xp + cond_ym + cond_yp + cond_zm + cond_zp)
            
            ss = cond_xm * p[i-1, j, k] + cond_xp * p[i+1, j, k] +
                 cond_ym * p[i, j-1, k] + cond_yp * p[i, j+1, k] +
                 cond_zm * p[i, j, k-1] + cond_zp * p[i, j, k+1]
            
            # Ap = A * p = -(ss - dd * p) = dd * p - ss
            Ap[i, j, k] = (dd * p[i, j, k] - ss) * m0
        end
    end
end


"""
    calc_residual_cg!(r, p, rhs, mask, grid, par)
    
Compute residual r = b - Ap and return ||r||
"""
function calc_residual_cg!(
    r::Array{Float64, 3},
    p::Array{Float64, 3},
    rhs::Array{Float64, 3},
    mask::Array{Float64, 3},
    grid::GridData,
    par::String
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
            
            dd = (1.0 - m0) + (cond_xm + cond_xp + cond_ym + cond_yp + cond_zm + cond_zp)
            
            ss = cond_xm * p[i-1, j, k] + cond_xp * p[i+1, j, k] +
                 cond_ym * p[i, j-1, k] + cond_yp * p[i, j+1, k] +
                 cond_zm * p[i, j, k-1] + cond_zp * p[i, j, k+1]
            
            b_val = rhs[i, j, k] * vol
            
            # r = b - Ap = b - (dd * p - ss) = b - dd * p + ss
            # Note: Our discretization is A*p = -∇²p, so r = b - A*p
            # Wait, need to carefully match sign convention with SOR.
            # In SOR: ss - b_val = dd * p  =>  (ss - b_val) / dd = p
            # So: A*p = dd*p - ss  (implicitly)
            # RHS convention: rhs = ∇·u* / dt, so b = rhs * vol
            # r = b - A*p = b - (dd*p - ss) = b + ss - dd*p
            r_val = (b_val + ss - dd * p[i, j, k]) * m0
            r[i, j, k] = r_val
            res_sq += r_val * r_val
        end
    end
    
    return sqrt(res_sq)
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
    solve_poisson_bicgstab!(buffers, grid, config, bc_set, par)
    
BiCGSTAB法による実装 (Placeholder)
"""
function solve_poisson_bicgstab!(
    buffers::CFDBuffers,
    grid::GridData,
    config::PoissonConfig,
    bc_set,
    par::String
)::Tuple{Bool, Int, Float64}
    # TODO: Implement BiCGSTAB solver
    # error("BiCGSTAB solver not implemented yet")
    println("WARNING: BiCGSTAB solver called but not implemented. Doing nothing.")
    return (false, 0, 1.0)
end

end # module PressureSolver
