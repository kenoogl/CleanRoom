module FractionalStep

using ..Common
using ..Grid
using ..Fields
using ..Convection
using ..Diffusion
using ..Turbulence
using ..PressureSolver
using ..BoundaryConditions

export fractional_step!, interpolate_to_faces!, compute_divergence!, correct_velocity!, compute_pseudo_velocity!

"""
    compute_pseudo_velocity!(buffers, grid, dt, par)

擬似速度 u* = u^n + dt * (Convection + Diffusion) を計算。
"""
function compute_pseudo_velocity!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    bc_set::BoundaryConditionSet,
    par::String
)
    # Initialize flux with 0
    # Actually add_convection/diffusion add to flux.
    # So we must clear flux first.
    fill!(buffers.flux_u, 0.0)
    fill!(buffers.flux_v, 0.0)
    fill!(buffers.flux_w, 0.0)
    
    # 1. Convection
    add_convection_flux!(buffers, grid, par)
    
    # 2. Diffusion
    add_diffusion_flux!(buffers, grid, bc_set, par)
    
    # 3. Update u*
    # u* = u + dt * flux
    # Only update fluid cells?
    # mask handling implied in fluxes (zero if wall).
    # so u* at wall = u_wall + dt*0 = u_wall (0).
    
    @inbounds @. buffers.u_star = buffers.u + dt * buffers.flux_u
    @inbounds @. buffers.v_star = buffers.v + dt * buffers.flux_v
    @inbounds @. buffers.w_star = buffers.w + dt * buffers.flux_w
end

"""
    interpolate_to_faces!(buffers, grid, par)

セル中心の擬似速度をセルフェイスへ内挿。
u_face_x[i] corresponds to u_{i-1/2} or u_{i+1/2}?
Grid convention: z_face[k] is between k-1 and k.
So u_face_x[i] should be at i-1/2 (left face of cell i).
"""
function interpolate_to_faces!(
    buffers::CFDBuffers,
    grid::GridData,
    par::String
)
    # Simple arithmetic mean: u_{i-1/2} = 0.5 * (u_{i-1} + u_i)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    # x-faces (range i=2 to mx ?)
    # u_face_x[i,j,k] = 0.5*(u*[i-1] + u*[i])
    @inbounds @. buffers.u_face_x[2:end, :, :] = 0.5 * (buffers.u_star[1:end-1, :, :] + buffers.u_star[2:end, :, :]) * 
                                                 buffers.mask[1:end-1, :, :] * buffers.mask[2:end, :, :]
    
    # y-faces
    @inbounds @. buffers.v_face_y[:, 2:end, :] = 0.5 * (buffers.v_star[:, 1:end-1, :] + buffers.v_star[:, 2:end, :]) * 
                                                 buffers.mask[:, 1:end-1, :] * buffers.mask[:, 2:end, :]
    
    # z-faces
    @inbounds @. buffers.w_face_z[:, :, 2:end] = 0.5 * (buffers.w_star[:, :, 1:end-1] + buffers.w_star[:, :, 2:end]) * 
                                                 buffers.mask[:, :, 1:end-1] * buffers.mask[:, :, 2:end]
    
    # Note: Boundary faces (ex i=3 start of real domain).
    # u_face_x[3] = 0.5*(u[2] + u[3]). u[2] is ghost.
    # Ghost values should be populated before this step?
    # Yes, BCs applied at end of previous step.
end

"""
    compute_divergence!(buffers, grid, dt, par)

発散計算とポアソン右辺生成。
rhs = (1/dt) * ∇·u*
"""
function compute_divergence!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    par::String
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    inv_dt = 1.0 / dt
    
    @inbounds for k in 3:mz-2
        dz = grid.dz[k]
        for j in 3:my-2
            for i in 3:mx-2
                if buffers.mask[i, j, k] == 0.0
                    buffers.rhs[i, j, k] = 0.0
                    continue
                end
                
                # Div u* = (u*_{i+1/2} - u*_{i-1/2})/dx + ...
                # u_face_x[i] is left face (i-1/2)?
                # u_face_x[i+1] is right face (i+1/2).
                
                div = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx +
                      (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy +
                      (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
                
                buffers.rhs[i, j, k] = div * inv_dt
            end
        end
    end
end

"""
    correct_velocity!(buffers, grid, dt, par)

速度補正。
1. セルフェイス速度を圧力勾配で修正: u_face^{n+1} = u_face* - dt * ∇p
2. セルセンター速度を両側の圧力勾配の平均で修正: u^{n+1} = u* - dt * avg(∇p)
"""
function correct_velocity!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    par::String
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    p = buffers.p
    mask = buffers.mask
    
    # 1. セルフェイス速度の修正
    # X方向フェイス (i-1/2位置)
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-1
        # 両側のセルが流体の場合のみ修正
        m_left = mask[i-1, j, k]
        m_right = mask[i, j, k]
        if m_left > 0 && m_right > 0
            dp_dx = (p[i, j, k] - p[i-1, j, k]) / dx
            buffers.u_face_x[i, j, k] -= dt * dp_dx
        else
            # 壁面の場合はゼロ
            buffers.u_face_x[i, j, k] = 0.0
        end
    end
    
    # Y方向フェイス (j-1/2位置)
    @inbounds for k in 3:mz-2, j in 3:my-1, i in 3:mx-2
        m_left = mask[i, j-1, k]
        m_right = mask[i, j, k]
        if m_left > 0 && m_right > 0
            dp_dy = (p[i, j, k] - p[i, j-1, k]) / dy
            buffers.v_face_y[i, j, k] -= dt * dp_dy
        else
            buffers.v_face_y[i, j, k] = 0.0
        end
    end
    
    # Z方向フェイス (k-1/2位置)
    @inbounds for k in 3:mz-1, j in 3:my-2, i in 3:mx-2
        m_left = mask[i, j, k-1]
        m_right = mask[i, j, k]
        if m_left > 0 && m_right > 0
            dz = grid.z_face[k] - grid.z_face[k-1]
            dp_dz = (p[i, j, k] - p[i, j, k-1]) / dz
            buffers.w_face_z[i, j, k] -= dt * dp_dz
        else
            buffers.w_face_z[i, j, k] = 0.0
        end
    end
    
    # 2. セルセンター速度の修正（両側フェイスの圧力勾配の平均）
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        if mask[i, j, k] == 0.0
            # 壁面セル
            buffers.u[i, j, k] = 0.0
            buffers.v[i, j, k] = 0.0
            buffers.w[i, j, k] = 0.0
            continue
        end
        
        # X方向: 左右フェイスの圧力勾配を平均
        m_left = mask[i-1, j, k]
        m_right = mask[i+1, j, k]
        
        # 左フェイス (i-1/2) の圧力勾配（マスクでNeumann条件）
        dp_dx_left = m_left * (p[i, j, k] - p[i-1, j, k]) / dx
        # 右フェイス (i+1/2) の圧力勾配
        dp_dx_right = m_right * (p[i+1, j, k] - p[i, j, k]) / dx
        
        # 平均勾配で修正
        dp_dx_avg = 0.5 * (dp_dx_left + dp_dx_right)
        buffers.u[i, j, k] = buffers.u_star[i, j, k] - dt * dp_dx_avg
        
        # Y方向
        m_front = mask[i, j-1, k]
        m_back = mask[i, j+1, k]
        
        dp_dy_front = m_front * (p[i, j, k] - p[i, j-1, k]) / dy
        dp_dy_back = m_back * (p[i, j+1, k] - p[i, j, k]) / dy
        
        dp_dy_avg = 0.5 * (dp_dy_front + dp_dy_back)
        buffers.v[i, j, k] = buffers.v_star[i, j, k] - dt * dp_dy_avg
        
        # Z方向
        m_bottom = mask[i, j, k-1]
        m_top = mask[i, j, k+1]
        
        dz_bottom = grid.z_face[k] - grid.z_face[k-1]
        dz_top = grid.z_face[k+1] - grid.z_face[k]
        
        dp_dz_bottom = m_bottom * (p[i, j, k] - p[i, j, k-1]) / dz_bottom
        dp_dz_top = m_top * (p[i, j, k+1] - p[i, j, k]) / dz_top
        
        dp_dz_avg = 0.5 * (dp_dz_bottom + dp_dz_top)
        buffers.w[i, j, k] = buffers.w_star[i, j, k] - dt * dp_dz_avg
    end
end

"""
    fractional_step!(buffers, grid, dt, bc_set, poisson_config, Cs, par)

統合ステップ関数。
"""
function fractional_step!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    bc_set::BoundaryConditionSet,
    poisson_config::PoissonConfig,
    Cs::Float64,
    nu::Float64,
    par::String
)::Tuple{Int, Float64}

    # 1. Update Turbulence (nu_t)
    compute_turbulent_viscosity!(buffers.nu_t, buffers, grid, Cs, nu, par)
    
    # 2. Pseudo Velocity
    compute_pseudo_velocity!(buffers, grid, dt, bc_set, par)
    
    # 2.5 Apply BCs to Pseudo Velocity
    apply_velocity_bcs!(buffers.u_star, buffers.v_star, buffers.w_star, grid, buffers.mask, bc_set, dt)

    # 3. Interpolate to Faces
    interpolate_to_faces!(buffers, grid, par)
    
    # 4. Divergence / RHS
    compute_divergence!(buffers, grid, dt, par)
    
    # 5. Poisson Solve
    converged, iter, res = solve_poisson!(buffers, grid, poisson_config, bc_set, par)
    
    # 6. Correct Velocity (Face and Center)
    correct_velocity!(buffers, grid, dt, par)
    
    # 7. Apply BCs to Face Velocities
    # セルフェイス速度への境界条件適用（周期境界、対称境界など）
    apply_face_velocity_bcs!(buffers.u_face_x, buffers.v_face_y, buffers.w_face_z, grid, buffers.mask, bc_set, dt)
    
    # 8. Apply BCs to Center Velocities
    apply_boundary_conditions!(buffers, grid, bc_set, dt, par)
    
    # 9. Update Time Average (Should be called by Main? Or here?)
    # Design flow says "Update Time Average" is last step of FracStep flow.
    # But usually caller manages averaging start time.
    # I'll leave it to Main because it depends on Time.
    
    return (iter, res)
end

end # module FractionalStep
