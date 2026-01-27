module Turbulence

using ..Common
using ..Grid
using ..Fields

export compute_turbulent_viscosity!

"""
    compute_turbulent_viscosity!(buffers, grid, Cs, nu_lam, par)

標準Smagorinskyモデルによる乱流粘性係数の計算。
νt = (Cs·Δ)² · |S|
todo compactなフィルタ幅の計算
"""
function compute_turbulent_viscosity!(
    buffers::CFDBuffers,
    grid::GridData,
    Cs::Float64,
    nu_lam::Float64,
    par::String
)
    # フィルタ幅 Δ = (Δx * Δy * Δz)^(1/3)
    # 各セルで計算（非等間隔格子のため）
    
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    # 歪み速度テンソル S_ij = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)
    # |S| = sqrt(2 * Σ S_ij S_ij)
    
    # ループ (i, j, k) 実領域 + ゴースト考慮？
    # 乱流粘性は全域で計算しておく（BC適用にも便利）
    
    # Parallel execution if par == "thread" (handled by FLoops in future)
    # For now straightforward loop
    
    dx = grid.dx
    dy = grid.dy
    
    @inbounds for k in 2:mz-1, j in 2:my-1, i in 2:mx-1
        dz = grid.dz[k]
        
        # 速度勾配を中心差分で
        dudx = (buffers.u[i+1, j, k] - buffers.u[i-1, j, k]) / (2*dx)
        dudy = (buffers.u[i, j+1, k] - buffers.u[i, j-1, k]) / (2*dy)
        dudz = (buffers.u[i, j, k+1] - buffers.u[i, j, k-1]) / (2*dz)
        
        dvdx = (buffers.v[i+1, j, k] - buffers.v[i-1, j, k]) / (2*dx)
        dvdy = (buffers.v[i, j+1, k] - buffers.v[i, j-1, k]) / (2*dy)
        dvdz = (buffers.v[i, j, k+1] - buffers.v[i, j, k-1]) / (2*dz)
        
        dwdx = (buffers.w[i+1, j, k] - buffers.w[i-1, j, k]) / (2*dx)
        dwdy = (buffers.w[i, j+1, k] - buffers.w[i, j-1, k]) / (2*dy)
        dwdz = (buffers.w[i, j, k+1] - buffers.w[i, j, k-1]) / (2*dz)
        
        S11 = dudx
        S22 = dvdy
        S33 = dwdz
        S12 = 0.5 * (dudy + dvdx)
        S13 = 0.5 * (dudz + dwdx)
        S23 = 0.5 * (dvdz + dwdy)
        
        # SijSij = S11^2 + S22^2 + S33^2 + 2*(S12^2 + S13^2 + 23^2)
        S_mag_sq = S11^2 + S22^2 + S33^2 + 2.0 * (S12^2 + S13^2 + S23^2)
        S_mag = sqrt(2.0 * S_mag_sq)
        
        delta = (dx * dy * dz)^(1.0/3.0)
        
        nu_t_val = (Cs * delta)^2 * S_mag
        buffers.nu_t[i, j, k] = nu_t_val
        buffers.nu_eff[i, j, k] = nu_lam + nu_t_val
    end
end

end # module Turbulence
