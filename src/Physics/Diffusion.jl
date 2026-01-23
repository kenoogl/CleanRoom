module Diffusion

using ..Common
using ..Grid
using ..Fields
using ..BoundaryConditions: BoundaryConditionSet, Wall, SlidingWall, Inlet

export add_diffusion_flux!

"""
    harmonic_mean(a, b)

調和平均を計算。
"""
@inline function harmonic_mean(a::Float64, b::Float64)
    return (a * b == 0.0) ? 0.0 : 2.0 * a * b / (a + b)
end

"""
    add_diffusion_flux!(buffers, grid, bc_set, par)

拡散フラックスを計算し、buffers.fluxに加算する。
2次精度中心差分、界面粘性は調和平均。
Wall/SlidingWallに対しては壁面せん断補正を適用する。
"""
function add_diffusion_flux!(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    par::String
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    dx2, dy2 = dx^2, dy^2
    
    # 1. 内部セル（3:mx-2...）の基本ループ（Neumann条件考慮済み）
    @inbounds for k in 3:mz-2
        dz = grid.dz[k] 
        dz2 = dz^2
        dz_m = grid.z_center[k] - grid.z_center[k-1] 
        dz_p = grid.z_center[k+1] - grid.z_center[k]
        
        for j in 3:my-2
            for i in 3:mx-2
                # --- セル中心の値 ---
                u_c = buffers.u[i, j, k]
                v_c = buffers.v[i, j, k]
                w_c = buffers.w[i, j, k]
                nu_c = buffers.nu_eff[i, j, k]
                
                # --- u-component ---
                # x-direction
                nu_w = harmonic_mean(buffers.nu_eff[i-1, j, k], nu_c)
                nu_e = harmonic_mean(nu_c, buffers.nu_eff[i+1, j, k])
                grad_w = (u_c - buffers.u[i-1, j, k]) / dx
                grad_e = (buffers.u[i+1, j, k] - u_c) / dx
                if buffers.mask[i-1, j, k] == 0.0; grad_w = 0.0; end
                if buffers.mask[i+1, j, k] == 0.0; grad_e = 0.0; end
                diff_u_x = (nu_e * grad_e - nu_w * grad_w) / dx
                
                # y-direction
                nu_s = harmonic_mean(buffers.nu_eff[i, j-1, k], nu_c)
                nu_n = harmonic_mean(nu_c, buffers.nu_eff[i, j+1, k])
                grad_s = (u_c - buffers.u[i, j-1, k]) / dy
                grad_n = (buffers.u[i, j+1, k] - u_c) / dy
                if buffers.mask[i, j-1, k] == 0.0; grad_s = 0.0; end
                if buffers.mask[i, j+1, k] == 0.0; grad_n = 0.0; end
                diff_u_y = (nu_n * grad_n - nu_s * grad_s) / dy
                
                # z-direction
                nu_b = harmonic_mean(buffers.nu_eff[i, j, k-1], nu_c)
                nu_t = harmonic_mean(nu_c, buffers.nu_eff[i, j, k+1])
                grad_b = (u_c - buffers.u[i, j, k-1]) / dz_m
                grad_t = (buffers.u[i, j, k+1] - u_c) / dz_p
                if buffers.mask[i, j, k-1] == 0.0; grad_b = 0.0; end
                if buffers.mask[i, j, k+1] == 0.0; grad_t = 0.0; end
                diff_u_z = (nu_t * grad_t - nu_b * grad_b) / dz
                
                buffers.flux_u[i, j, k] += (diff_u_x + diff_u_y + diff_u_z)
                
                # --- v-component ---
                # x
                grad_w = (v_c - buffers.v[i-1, j, k]) / dx
                grad_e = (buffers.v[i+1, j, k] - v_c) / dx
                if buffers.mask[i-1, j, k] == 0.0; grad_w = 0.0; end
                if buffers.mask[i+1, j, k] == 0.0; grad_e = 0.0; end
                diff_v_x = (nu_e * grad_e - nu_w * grad_w) / dx
                
                # y
                grad_s = (v_c - buffers.v[i, j-1, k]) / dy
                grad_n = (buffers.v[i, j+1, k] - v_c) / dy
                if buffers.mask[i, j-1, k] == 0.0; grad_s = 0.0; end
                if buffers.mask[i, j+1, k] == 0.0; grad_n = 0.0; end
                diff_v_y = (nu_n * grad_n - nu_s * grad_s) / dy
                
                # z
                grad_b = (v_c - buffers.v[i, j, k-1]) / dz_m
                grad_t = (buffers.v[i, j, k+1] - v_c) / dz_p
                if buffers.mask[i, j, k-1] == 0.0; grad_b = 0.0; end
                if buffers.mask[i, j, k+1] == 0.0; grad_t = 0.0; end
                diff_v_z = (nu_t * grad_t - nu_b * grad_b) / dz
                
                buffers.flux_v[i, j, k] += (diff_v_x + diff_v_y + diff_v_z)
                
                # --- w-component ---
                # x
                grad_w = (w_c - buffers.w[i-1, j, k]) / dx
                grad_e = (buffers.w[i+1, j, k] - w_c) / dx
                if buffers.mask[i-1, j, k] == 0.0; grad_w = 0.0; end
                if buffers.mask[i+1, j, k] == 0.0; grad_e = 0.0; end
                diff_w_x = (nu_e * grad_e - nu_w * grad_w) / dx
                
                # y
                grad_s = (w_c - buffers.w[i, j-1, k]) / dy
                grad_n = (buffers.w[i, j+1, k] - w_c) / dy
                if buffers.mask[i, j-1, k] == 0.0; grad_s = 0.0; end
                if buffers.mask[i, j+1, k] == 0.0; grad_n = 0.0; end
                diff_w_y = (nu_n * grad_n - nu_s * grad_s) / dy
                
                # z
                grad_b = (w_c - buffers.w[i, j, k-1]) / dz_m
                grad_t = (buffers.w[i, j, k+1] - w_c) / dz_p
                if buffers.mask[i, j, k-1] == 0.0; grad_b = 0.0; end
                if buffers.mask[i, j, k+1] == 0.0; grad_t = 0.0; end
                diff_w_z = (nu_t * grad_t - nu_b * grad_b) / dz
                
                buffers.flux_w[i, j, k] += (diff_w_x + diff_w_y + diff_w_z)
            end
        end
    end

    # 2. 壁面せん断補正（Wall / SlidingWall）
    # ループ内で Neumann (grad=0) になっている箇所に、2*nu_eff*(Uw - u)/dn^2 を加算
    
    # 外部境界の各面をチェック
    faces = [(:x_min, bc_set.x_min), (:x_max, bc_set.x_max), 
             (:y_min, bc_set.y_min), (:y_max, bc_set.y_max), 
             (:z_min, bc_set.z_min), (:z_max, bc_set.z_max)]
    
    for (f_sym, bc) in faces
        if bc.velocity_type == Wall || bc.velocity_type == SlidingWall || bc.velocity_type == Inlet
            uw, vw, ww = bc.velocity_value
            nu_eff = buffers.nu_eff
            
            if f_sym == :x_min
                i = 3
                @inbounds for k in 3:mz-2, j in 3:my-2
                    buffers.flux_u[i,j,k] += 2.0 * nu_eff[i,j,k] * (uw - buffers.u[i,j,k]) / dx2
                    buffers.flux_v[i,j,k] += 2.0 * nu_eff[i,j,k] * (vw - buffers.v[i,j,k]) / dx2
                    buffers.flux_w[i,j,k] += 2.0 * nu_eff[i,j,k] * (ww - buffers.w[i,j,k]) / dx2
                end
            elseif f_sym == :x_max
                i = mx-2
                @inbounds for k in 3:mz-2, j in 3:my-2
                    buffers.flux_u[i,j,k] += 2.0 * nu_eff[i,j,k] * (uw - buffers.u[i,j,k]) / dx2
                    buffers.flux_v[i,j,k] += 2.0 * nu_eff[i,j,k] * (vw - buffers.v[i,j,k]) / dx2
                    buffers.flux_w[i,j,k] += 2.0 * nu_eff[i,j,k] * (ww - buffers.w[i,j,k]) / dx2
                end
            elseif f_sym == :y_min
                j = 3
                @inbounds for k in 3:mz-2, i in 3:mx-2
                    buffers.flux_u[i,j,k] += 2.0 * nu_eff[i,j,k] * (uw - buffers.u[i,j,k]) / dy2
                    buffers.flux_v[i,j,k] += 2.0 * nu_eff[i,j,k] * (vw - buffers.v[i,j,k]) / dy2
                    buffers.flux_w[i,j,k] += 2.0 * nu_eff[i,j,k] * (ww - buffers.w[i,j,k]) / dy2
                end
            elseif f_sym == :y_max
                j = my-2
                @inbounds for k in 3:mz-2, i in 3:mx-2
                    buffers.flux_u[i,j,k] += 2.0 * nu_eff[i,j,k] * (uw - buffers.u[i,j,k]) / dy2
                    buffers.flux_v[i,j,k] += 2.0 * nu_eff[i,j,k] * (vw - buffers.v[i,j,k]) / dy2
                    buffers.flux_w[i,j,k] += 2.0 * nu_eff[i,j,k] * (ww - buffers.w[i,j,k]) / dy2
                end
            elseif f_sym == :z_min
                k = 3
                dz2 = grid.dz[k]^2
                @inbounds for j in 3:my-2, i in 3:mx-2
                    buffers.flux_u[i,j,k] += 2.0 * nu_eff[i,j,k] * (uw - buffers.u[i,j,k]) / dz2
                    buffers.flux_v[i,j,k] += 2.0 * nu_eff[i,j,k] * (vw - buffers.v[i,j,k]) / dz2
                    buffers.flux_w[i,j,k] += 2.0 * nu_eff[i,j,k] * (ww - buffers.w[i,j,k]) / dz2
                end
            elseif f_sym == :z_max
                k = mz-2
                dz2 = grid.dz[k]^2
                @inbounds for j in 3:my-2, i in 3:mx-2
                    buffers.flux_u[i,j,k] += 2.0 * nu_eff[i,j,k] * (uw - buffers.u[i,j,k]) / dz2
                    buffers.flux_v[i,j,k] += 2.0 * nu_eff[i,j,k] * (vw - buffers.v[i,j,k]) / dz2
                    buffers.flux_w[i,j,k] += 2.0 * nu_eff[i,j,k] * (ww - buffers.w[i,j,k]) / dz2
                end
            end
        end
    end
end

end # module Diffusion
