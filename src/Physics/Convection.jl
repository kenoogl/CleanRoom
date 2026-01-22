module Convection

using ..Common
using ..Grid
using ..Fields

export add_convection_flux!

"""
    weno3_reconstruct_left(u_im1, u_i, u_ip1, dx_im1, dx_i, dx_ip1, epsilon)

WENO3再構成（左側 u^{-}_{i+1/2}）。
"""
@inline function weno3_reconstruct_left(
    u_im1::Float64, u_i::Float64, u_ip1::Float64,
    dx_im1::Float64, dx_i::Float64, dx_ip1::Float64,
    epsilon::Float64
)
    # 簡易実装: 等間隔係数を使用
    p0 = -0.5*u_im1 + 1.5*u_i
    p1 = 0.5*u_i + 0.5*u_ip1
    
    beta0 = (u_i - u_im1)^2
    beta1 = (u_ip1 - u_i)^2
    
    alpha0 = (1.0/3.0) / (epsilon + beta0)^2
    alpha1 = (2.0/3.0) / (epsilon + beta1)^2
    
    omega0 = alpha0 / (alpha0 + alpha1)
    omega1 = alpha1 / (alpha0 + alpha1)
    
    return omega0 * p0 + omega1 * p1
end

"""
    add_convection_flux!(buffers, grid, par)

Lax-Friedrichs flux splitting + WENO3 reconstruction.
"""
function add_convection_flux!(
    buffers::CFDBuffers,
    grid::GridData,
    par::String
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    # --- Advection of U ---
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 2:mx-2
        # Flux F = u*u in x
        u_i, u_ip1 = buffers.u[i, j, k], buffers.u[i+1, j, k]
        alpha = max(abs(u_i), abs(u_ip1))
        
        fp(u) = 0.5 * (u*u + alpha*u)
        fm(u) = 0.5 * (u*u - alpha*u)
        
        val_p = weno3_reconstruct_left(fp(buffers.u[i-1, j, k]), fp(u_i), fp(u_ip1), grid.dx, grid.dx, grid.dx, 1e-6)
        val_m = weno3_reconstruct_left(fm(buffers.u[i+2, j, k]), fm(u_ip1), fm(u_i), grid.dx, grid.dx, grid.dx, 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i+1, j, k])
        buffers.flux_u[i, j, k] -= flx / grid.dx
        buffers.flux_u[i+1, j, k] += flx / grid.dx
    end

    @inbounds for k in 3:mz-2, j in 2:my-2, i in 3:mx-2
        # Flux F = v*u in y
        v_j, v_jp1 = buffers.v[i, j, k], buffers.v[i, j+1, k]
        alpha = max(abs(v_j), abs(v_jp1))
        
        fp_y(u, v) = 0.5 * (v*u + alpha*u)
        fm_y(u, v) = 0.5 * (v*u - alpha*u)
        
        val_p = weno3_reconstruct_left(fp_y(buffers.u[i, j-1, k], buffers.v[i, j-1, k]), fp_y(buffers.u[i, j, k], v_j), fp_y(buffers.u[i, j+1, k], v_jp1), grid.dy, grid.dy, grid.dy, 1e-6)
        val_m = weno3_reconstruct_left(fm_y(buffers.u[i, j+2, k], buffers.v[i, j+2, k]), fm_y(buffers.u[i, j+1, k], v_jp1), fm_y(buffers.u[i, j, k], v_j), grid.dy, grid.dy, grid.dy, 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i, j+1, k])
        buffers.flux_u[i, j, k] -= flx / grid.dy
        buffers.flux_u[i, j+1, k] += flx / grid.dy
    end

    @inbounds for k in 2:mz-2, j in 3:my-2, i in 3:mx-2
        # Flux F = w*u in z
        w_k, w_kp1 = buffers.w[i, j, k], buffers.w[i, j, k+1]
        alpha = max(abs(w_k), abs(w_kp1))
        
        fp_z(u, w) = 0.5 * (w*u + alpha*u)
        fm_z(u, w) = 0.5 * (w*u - alpha*u)
        
        val_p = weno3_reconstruct_left(fp_z(buffers.u[i, j, k-1], buffers.w[i, j, k-1]), fp_z(buffers.u[i, j, k], w_k), fp_z(buffers.u[i, j, k+1], w_kp1), grid.dz[k-1], grid.dz[k], grid.dz[k+1], 1e-6)
        val_m = weno3_reconstruct_left(fm_z(buffers.u[i, j, k+2], buffers.w[i, j, k+2]), fm_z(buffers.u[i, j, k+1], w_kp1), fm_z(buffers.u[i, j, k], w_k), grid.dz[k+2], grid.dz[k+1], grid.dz[k], 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i, j, k+1])
        buffers.flux_u[i, j, k] -= flx / grid.dz[k]
        buffers.flux_u[i, j, k+1] += flx / grid.dz[k+1]
    end

    # --- Advection of V ---
    # similar structure for flux_v

    @inbounds for k in 3:mz-2, j in 3:my-2, i in 2:mx-2
        # F = u*v in x
        u_i, u_ip1 = buffers.u[i, j, k], buffers.u[i+1, j, k]
        alpha = max(abs(u_i), abs(u_ip1))
        fp_x(v, u) = 0.5 * (u*v + alpha*v)
        fm_x(v, u) = 0.5 * (u*v - alpha*v)
        
        val_p = weno3_reconstruct_left(fp_x(buffers.v[i-1, j, k], buffers.u[i-1, j, k]), fp_x(buffers.v[i, j, k], u_i), fp_x(buffers.v[i+1, j, k], u_ip1), grid.dx, grid.dx, grid.dx, 1e-6)
        val_m = weno3_reconstruct_left(fm_x(buffers.v[i+2, j, k], buffers.u[i+2, j, k]), fm_x(buffers.v[i+1, j, k], u_ip1), fm_x(buffers.v[i, j, k], u_i), grid.dx, grid.dx, grid.dx, 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i+1, j, k])
        buffers.flux_v[i, j, k] -= flx / grid.dx
        buffers.flux_v[i+1, j, k] += flx / grid.dx
    end
    
    @inbounds for k in 3:mz-2, j in 2:my-2, i in 3:mx-2
        # F = v*v in y
        v_j, v_jp1 = buffers.v[i, j, k], buffers.v[i, j+1, k]
        alpha = max(abs(v_j), abs(v_jp1))
        fp(v) = 0.5 * (v*v + alpha*v)
        fm(v) = 0.5 * (v*v - alpha*v)
        
        val_p = weno3_reconstruct_left(fp(buffers.v[i, j-1, k]), fp(v_j), fp(v_jp1), grid.dy, grid.dy, grid.dy, 1e-6)
        val_m = weno3_reconstruct_left(fm(buffers.v[i, j+2, k]), fm(v_jp1), fm(v_j), grid.dy, grid.dy, grid.dy, 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i, j+1, k])
        buffers.flux_v[i, j, k] -= flx / grid.dy
        buffers.flux_v[i, j+1, k] += flx / grid.dy
    end

    @inbounds for k in 2:mz-2, j in 3:my-2, i in 3:mx-2
        # F = w*v in z
        w_k, w_kp1 = buffers.w[i, j, k], buffers.w[i, j, k+1]
        alpha = max(abs(w_k), abs(w_kp1))
        
        fp_z(v, w) = 0.5 * (w*v + alpha*v)
        fm_z(v, w) = 0.5 * (w*v - alpha*v)
        
        val_p = weno3_reconstruct_left(fp_z(buffers.v[i, j, k-1], buffers.w[i, j, k-1]), fp_z(buffers.v[i, j, k], w_k), fp_z(buffers.v[i, j, k+1], w_kp1), grid.dz[k-1], grid.dz[k], grid.dz[k+1], 1e-6)
        val_m = weno3_reconstruct_left(fm_z(buffers.v[i, j, k+2], buffers.w[i, j, k+2]), fm_z(buffers.v[i, j, k+1], w_kp1), fm_z(buffers.v[i, j, k], w_k), grid.dz[k+2], grid.dz[k+1], grid.dz[k], 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i, j, k+1])
        buffers.flux_v[i, j, k] -= flx / grid.dz[k]
        buffers.flux_v[i, j, k+1] += flx / grid.dz[k+1]
    end

    # --- Advection of W ---
    # similar structure for flux_w

       @inbounds for k in 3:mz-2, j in 3:my-2, i in 2:mx-2
        # F = u*w in x
        u_i, u_ip1 = buffers.u[i, j, k], buffers.u[i+1, j, k]
        alpha = max(abs(u_i), abs(u_ip1))
        
        fp_x(w, u) = 0.5 * (u*w + alpha*w)
        fm_x(w, u) = 0.5 * (u*w - alpha*w)
        
        val_p = weno3_reconstruct_left(fp_x(buffers.w[i-1, j, k], buffers.u[i-1, j, k]), fp_x(buffers.w[i, j, k], u_i), fp_x(buffers.w[i+1, j, k], u_ip1), grid.dx, grid.dx, grid.dx, 1e-6)
        val_m = weno3_reconstruct_left(fm_x(buffers.w[i+2, j, k], buffers.u[i+2, j, k]), fm_x(buffers.w[i+1, j, k], u_ip1), fm_x(buffers.w[i, j, k], u_i), grid.dx, grid.dx, grid.dx, 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i+1, j, k])
        buffers.flux_w[i, j, k] -= flx / grid.dx
        buffers.flux_w[i+1, j, k] += flx / grid.dx
    end

    @inbounds for k in 3:mz-2, j in 2:my-2, i in 3:mx-2
        # F = v*w in y
        v_j, v_jp1 = buffers.v[i, j, k], buffers.v[i, j+1, k]
        alpha = max(abs(v_j), abs(v_jp1))
        
        fp_y(w, v) = 0.5 * (v*w + alpha*w)
        fm_y(w, v) = 0.5 * (v*w - alpha*w)
        
        val_p = weno3_reconstruct_left(fp_y(buffers.w[i, j-1, k], buffers.v[i, j-1, k]), fp_y(buffers.w[i, j, k], v_j), fp_y(buffers.w[i, j+1, k], v_jp1), grid.dy, grid.dy, grid.dy, 1e-6)
        val_m = weno3_reconstruct_left(fm_y(buffers.w[i, j+2, k], buffers.v[i, j+2, k]), fm_y(buffers.w[i, j+1, k], v_jp1), fm_y(buffers.w[i, j, k], v_j), grid.dy, grid.dy, grid.dy, 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i, j+1, k])
        buffers.flux_w[i, j, k] -= flx / grid.dy
        buffers.flux_w[i, j+1, k] += flx / grid.dy
    end

    @inbounds for k in 2:mz-2, j in 3:my-2, i in 3:mx-2
        # F = w*w in z
        w_k, w_kp1 = buffers.w[i, j, k], buffers.w[i, j, k+1]
        alpha = max(abs(w_k), abs(w_kp1))
        
        fp(w) = 0.5 * (w*w + alpha*w)
        fm(w) = 0.5 * (w*w - alpha*w)
        
        val_p = weno3_reconstruct_left(fp(buffers.w[i, j, k-1]), fp(w_k), fp(w_kp1), grid.dz[k-1], grid.dz[k], grid.dz[k+1], 1e-6)
        val_m = weno3_reconstruct_left(fm(buffers.w[i, j, k+2]), fm(w_kp1), fm(w_k), grid.dz[k+2], grid.dz[k+1], grid.dz[k], 1e-6)
        
        flx = (val_p + val_m) * (buffers.mask[i, j, k] * buffers.mask[i, j, k+1])
        buffers.flux_w[i, j, k] -= flx / grid.dz[k]
        buffers.flux_w[i, j, k+1] += flx / grid.dz[k+1]
    end
end

end # module Convection
