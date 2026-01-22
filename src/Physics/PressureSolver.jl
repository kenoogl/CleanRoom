module PressureSolver

using ..Common
using ..Fields
using ..Grid
using ..BoundaryConditions
using LinearAlgebra # for norm

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
    idx2 = 1.0 / dx^2
    idy2 = 1.0 / dy^2
    
    while iter < config.max_iter
        iter += 1
        max_res = 0.0
        
        # SOR Loop (Red/Black)
        for color in 0:1
            @inbounds for k in 3:mz-2
                dZ_p = grid.z_center[k+1] - grid.z_center[k]
                dZ_m = grid.z_center[k] - grid.z_center[k-1]
                dz = grid.dz[k]
                
                c_t = 1.0 / (dz * dZ_p)
                c_b = 1.0 / (dz * dZ_m)
                
                for j in 3:my-2
                    for i in 3:mx-2
                        if mask[i, j, k] == 0.0
                            continue
                        end
                        
                        if (i + j + k) % 2 != color
                            continue
                        end
                        
                        denom = 0.0
                        sum_neighbors = 0.0
                        
                        # X-neighbors
                        if mask[i-1, j, k] > 0
                            sum_neighbors += p[i-1, j, k] * idx2
                            denom += idx2
                        end
                        if mask[i+1, j, k] > 0
                            sum_neighbors += p[i+1, j, k] * idx2
                            denom += idx2
                        end
                        
                        # Y-neighbors
                        if mask[i, j-1, k] > 0
                            sum_neighbors += p[i, j-1, k] * idy2
                            denom += idy2
                        end
                        if mask[i, j+1, k] > 0
                            sum_neighbors += p[i, j+1, k] * idy2
                            denom += idy2
                        end
                        
                        # Z-neighbors
                        if mask[i, j, k-1] > 0
                            sum_neighbors += p[i, j, k-1] * c_b
                            denom += c_b
                        end
                        if mask[i, j, k+1] > 0
                            sum_neighbors += p[i, j, k+1] * c_t
                            denom += c_t
                        end
                        
                        if denom == 0.0
                             continue
                        end
                        
                        # SOR update
                        p_star = (sum_neighbors - rhs[i, j, k]) / denom
                        p_new = (1.0 - config.omega) * p[i, j, k] + config.omega * p_star
                        
                        # Residual
                        res = abs(p_new - p[i, j, k])
                        if res > max_res
                            max_res = res
                        end
                        
                        p[i, j, k] = p_new
                    end
                end
            end
        end
        
        # Apply periodic BC
        if !isnothing(bc_set)
            apply_periodic_pressure!(p, grid, bc_set)
        end
        
        residual = max_res
        if residual < config.tol
            converged = true
            break
        end
    end
    
    return (converged, iter, residual)
end

"""
    solve_poisson_cg!(buffers, grid, config, bc_set, par)
    
Conjugate Gradient法による実装 (Placeholder)
"""
function solve_poisson_cg!(
    buffers::CFDBuffers,
    grid::GridData,
    config::PoissonConfig,
    bc_set,
    par::String
)::Tuple{Bool, Int, Float64}
    # TODO: Implement CG solver
    # error("CG solver not implemented yet")
    println("WARNING: CG solver called but not implemented. Doing nothing.")
    return (false, 0, 1.0)
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
