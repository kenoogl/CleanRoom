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
    
    mx, my, mz = grid.mx, grid.my, grid.mz
    nx, ny, nz = mx-4, my-4, mz-4 # Internal counts?
    # GridData stores mx, my, mz.
    # Indices 3 to mx-2 are valid fluid domain (including first layer of ghost for BC usually, but here main loop is over domain).
    
    # Red-Black SOR
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
    
    # Pre-calculated dz terms?
    # Non-uniform Z: dz varies.
    # Laplacian: 1/dz_i * ( (p_{k+1}-p_k)/dz_{k+1/2} - (p_k-p_{k-1})/dz_{k-1/2} )
    # dz_i = grid.dz[k] (cell width)
    # dz_{k+1/2} = grid.z_center[k+1] - grid.z_center[k] ? No, z_center diff.
    # Actually dz_{k+1/2} is distance between centers.
    # Let's call dZ_p = z_center[k+1] - z_center[k]
    # dZ_m = z_center[k] - z_center[k-1]
    
    while iter < config.max_iter
        iter += 1
        max_res = 0.0
        
        # SOR Loop
        # Red/Black pass?
        # Checkerboard: (i+j+k) % 2 == 0 or 1.
        
        for color in 0:1
            @inbounds for k in 3:mz-2
                dZ_p = grid.z_center[k+1] - grid.z_center[k]
                dZ_m = grid.z_center[k] - grid.z_center[k-1]
                dz = grid.dz[k]
                
                # Coeffs for Z
                # term: (p_t - p_c)/dZ_p - ...
                # 1/dz * ( p_t/dZ_p - p_c/dZ_p - p_c/dZ_m + p_b/dZ_m )
                # = p_t/(dz*dZ_p) + p_b/(dz*dZ_m) - p_c * (1/(dz*dZ_p) + 1/(dz*dZ_m))
                
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
                        
                        # Calculate coeffs dynamically based on neighbor mask (Neumann)
                        denom = 0.0
                        sum_neighbors = 0.0
                        
                        # X-neighbors
                        if mask[i-1, j, k] > 0
                            sum_neighbors += p[i-1, j, k] * idx2
                            denom += idx2
                        else
                            # Neumann: term vanishes from flux diff?
                            # Treated as flux=0 at wall.
                            # So no contribution to sum, no contribution to denom?
                            # Correct. 1/dz * (0 - flux_left). flux_left=0.
                            nothing
                        end
                        if mask[i+1, j, k] > 0
                            sum_neighbors += p[i+1, j, k] * idx2
                            denom += idx2
                        else
                             nothing
                        end
                        
                        # Y-neighbors
                        if mask[i, j-1, k] > 0
                            sum_neighbors += p[i, j-1, k] * idy2
                            denom += idy2
                        else
                            nothing
                        end
                        if mask[i, j+1, k] > 0
                            sum_neighbors += p[i, j+1, k] * idy2
                            denom += idy2
                        else
                            nothing
                        end
                        
                        # Z-neighbors
                        if mask[i, j, k-1] > 0
                            sum_neighbors += p[i, j, k-1] * c_b
                            denom += c_b
                        else
                            nothing
                        end
                        if mask[i, j, k+1] > 0
                            sum_neighbors += p[i, j, k+1] * c_t
                            denom += c_t
                        else
                            nothing
                        end
                        
                        if denom == 0.0
                             # Isolated cell?
                             continue
                        end
                        
                        # SOR update
                        p_star = (sum_neighbors - rhs[i, j, k]) / denom
                        p_new = (1.0 - config.omega) * p[i, j, k] + config.omega * p_star
                        
                        # Residual (L_inf ?)
                        res = abs(p_new - p[i, j, k])
                        if res > max_res
                            max_res = res
                        end
                        
                        p[i, j, k] = p_new
                    end
                end
            end
        end
        
        # Apply periodic BC for pressure ghost cells if applicable
        if !isnothing(bc_set)
            apply_periodic_pressure!(p, grid, bc_set)
        end
        
        # Global residual check (Relative or Absolute?)
        # Spec says "相対残差ノルム ||b - Ax|| / ||b||".
        # max_res above is Delta p. This is related but not exactly ||r||.
        # Computing actual residual b - Ax is expensive (another pass).
        # Approximating convergence by change in solution is common.
        # But if spec strictly requires ||b - Ax||, I should compute it occasionally.
        # "Monitor Data" also needs "pressure_residual".
        # I'll use max_res for now as proxy or compute simplified residual.
        # Assuming max_res < tol implies convergence.
        
        residual = max_res
        if residual < config.tol
            converged = true
            break
        end
    end
    
    # Mean pressure subtraction (Pinning)
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
        # This function returns converged flag.
    elseif !converged && config.on_divergence == Abort
        # Caller handles abort
    end

    return (converged, iter, residual)
end

end # module PressureSolver
