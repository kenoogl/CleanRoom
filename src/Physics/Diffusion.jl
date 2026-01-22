module Diffusion

using ..Common
using ..Grid
using ..Fields

export add_diffusion_flux!

"""
    harmonic_mean(a, b)

調和平均を計算。
"""
@inline function harmonic_mean(a::Float64, b::Float64)
    return (a * b == 0.0) ? 0.0 : 2.0 * a * b / (a + b)
end

"""
    add_diffusion_flux!(buffers, grid, par)

拡散フラックスを計算し、buffers.fluxに加算する。
2次精度中心差分、界面粘性は調和平均。
"""
function add_diffusion_flux!(
    buffers::CFDBuffers,
    grid::GridData,
    par::String
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    # u component diffusion
    # ∂/∂x(nu ∂u/∂x) + ∂/∂y(nu ∂u/∂y) + ∂/∂z(nu ∂u/∂z)
    
    # Loop over internal cells (3 to mx-2 etc if ghost=2)
    # Range: 3:mx-2 means real cells?
    # Spec: 5 point stencil (WENO) needs 2 ghosts.
    # Central diff needs 1 ghost.
    
    @inbounds for k in 3:mz-2
        dz = grid.dz[k] # Cell width
        # Distances between centers for gradients at faces
        
        # z-face k-0.5 (between k-1 and k)
        dz_m = grid.z_center[k] - grid.z_center[k-1] 
        # z-face k+0.5 (between k and k+1)
        dz_p = grid.z_center[k+1] - grid.z_center[k]
        
        for j in 3:my-2
            dy = grid.dy
            # Uniform grid: center dist = dy
            
            for i in 3:mx-2
                dx = grid.dx
                # Uniform grid: center dist = dx
                
                # --- u-component ---
                u_c = buffers.u[i, j, k]
                nu_c = buffers.nu_eff[i, j, k]
                
                # x-direction
                nu_w = harmonic_mean(buffers.nu_eff[i-1, j, k], nu_c)
                nu_e = harmonic_mean(nu_c, buffers.nu_eff[i+1, j, k])
                grad_w = (u_c - buffers.u[i-1, j, k]) / dx
                grad_e = (buffers.u[i+1, j, k] - u_c) / dx
                # mask check for wall. If mask=0 (solid), grad=0 (Neumann)
                # Apply to grad (or flux). 
                # Task 9.2: "壁面境界でのマスク関数による勾配ゼロ処理"
                # "u_ghost = u_i * mask_i + u_neighbor * (1 - mask_neighbor)" approach implies ghost value setting.
                # If we use explicit flux Calc:
                # If neighbor is solid, flux = 0? No, Neumann means ∂u/∂n = 0, so flux = 0 (since grad=0).
                # So if mask[i-1] == 0, grad_w = 0.
                
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

                diff_u_z = (nu_t * grad_t - nu_b * grad_b) / dz # divide by cell width dz
                
                buffers.flux_u[i, j, k] += (diff_u_x + diff_u_y + diff_u_z)
                
                # --- v-component ---
                v_c = buffers.v[i, j, k]
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
                w_c = buffers.w[i, j, k]
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
end

end # module Diffusion
