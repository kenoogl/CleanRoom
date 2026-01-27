module FractionalStep

using ..Common
using ..Grid
using ..Fields
using ..Convection
using ..Diffusion
using ..Turbulence
using ..PressureSolver
using ..BoundaryConditions

export fractional_step!, interpolate_to_faces!, compute_divergence!, correct_velocity!, compute_pseudo_velocity!, compute_conv_diff_flux!

"""
    compute_conv_diff_flux!(buffers, grid, bc_set, par)

対流+拡散のフラックスを計算（buffers.flux_*を更新）。
"""
function compute_conv_diff_flux!(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    par::String
)
    fill!(buffers.flux_u, 0.0)
    fill!(buffers.flux_v, 0.0)
    fill!(buffers.flux_w, 0.0)

    add_convection_flux!(buffers, grid, bc_set, par)
    add_diffusion_flux!(buffers, grid, bc_set, par)
end

"""
    compute_pseudo_velocity!(buffers, grid, dt, bc_set, par)

擬似速度 u* = u^n + dt * (Convection + Diffusion) を計算。
"""
function compute_pseudo_velocity!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    bc_set::BoundaryConditionSet,
    par::String
)
    compute_conv_diff_flux!(buffers, grid, bc_set, par)

    # u* = u + dt * flux
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
    apply_outflow_face_star!(buffers, grid, bc_set, dt, u_face_prev, v_face_prev, w_face_prev, reverse_flow_stabilization)

Outflow境界のセルフェイス擬似速度を、セルセンターとセルフェイス値から構成する。
u*_{i+1/2} = u^n_{i+1/2} - 2 Δt u_c (u^n_{i+1/2} - u^n_i) / Δx
"""
function apply_outflow_face_star!(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    dt::Float64,
    u_face_prev::Array{Float64, 3},
    v_face_prev::Array{Float64, 3},
    w_face_prev::Array{Float64, 3},
    reverse_flow_stabilization::Bool
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    u, v, w = buffers.u, buffers.v, buffers.w
    mask = buffers.mask

    uc = BoundaryConditions.compute_outflow_uc(u, v, w, grid, mask, bc_set)
    # clip_uc_p: for Max boundaries (outflow is u > 0). Clip if u < 0 (backflow).
    clip_uc_p(x) = (reverse_flow_stabilization && x < 0.0) ? 0.0 : x
    # clip_uc_n: for Min boundaries (outflow is u < 0). Clip if u > 0 (backflow). Return speed (positive).
    clip_uc_n(x) = (reverse_flow_stabilization && x > 0.0) ? 0.0 : -x

    if bc_set.x_max.velocity_type == Outflow
        i = mx - 2
        fi = mx - 1
        uc_face = clip_uc_p(uc.x_max)
        coeff = 2.0 * dt * uc_face / dx
        @inbounds for k in 3:mz-2, j in 3:my-2
            m0 = mask[i, j, k] * mask[fi, j, k]
            u_face_n = u_face_prev[fi, j, k]
            buffers.u_face_x[fi, j, k] = m0 * (u_face_n - coeff * (u_face_n - u[i, j, k]))
        end
    end
    if bc_set.x_min.velocity_type == Outflow
        i = 3
        fi = 3
        uc_face = clip_uc_n(uc.x_min)
        coeff = 2.0 * dt * uc_face / dx
        @inbounds for k in 3:mz-2, j in 3:my-2
            m0 = mask[i, j, k] * mask[fi-1, j, k]
            u_face_n = u_face_prev[fi, j, k]
            buffers.u_face_x[fi, j, k] = m0 * (u_face_n - coeff * (u_face_n - u[i, j, k]))
        end
    end

    if bc_set.y_max.velocity_type == Outflow
        j = my - 2
        fj = my - 1
        uc_face = clip_uc_p(uc.y_max)
        coeff = 2.0 * dt * uc_face / dy
        @inbounds for k in 3:mz-2, i in 3:mx-2
            m0 = mask[i, j, k] * mask[i, fj, k]
            v_face_n = v_face_prev[i, fj, k]
            buffers.v_face_y[i, fj, k] = m0 * (v_face_n - coeff * (v_face_n - v[i, j, k]))
        end
    end
    if bc_set.y_min.velocity_type == Outflow
        j = 3
        fj = 3
        uc_face = clip_uc_n(uc.y_min)
        coeff = 2.0 * dt * uc_face / dy
        @inbounds for k in 3:mz-2, i in 3:mx-2
            m0 = mask[i, j, k] * mask[i, fj-1, k]
            v_face_n = v_face_prev[i, fj, k]
            buffers.v_face_y[i, fj, k] = m0 * (v_face_n - coeff * (v_face_n - v[i, j, k]))
        end
    end

    if bc_set.z_max.velocity_type == Outflow
        k = mz - 2
        fk = mz - 1
        dz = grid.z_face[k+1] - grid.z_face[k]
        uc_face = clip_uc_p(uc.z_max)
        coeff = 2.0 * dt * uc_face / dz
        @inbounds for j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k] * mask[i, j, fk]
            w_face_n = w_face_prev[i, j, fk]
            buffers.w_face_z[i, j, fk] = m0 * (w_face_n - coeff * (w_face_n - w[i, j, k]))
        end
    end
    if bc_set.z_min.velocity_type == Outflow
        k = 3
        fk = 3
        dz = grid.z_face[k+1] - grid.z_face[k]
        uc_face = clip_uc_n(uc.z_min)
        coeff = 2.0 * dt * uc_face / dz
        @inbounds for j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k] * mask[i, j, fk-1]
            w_face_n = w_face_prev[i, j, fk]
            buffers.w_face_z[i, j, fk] = m0 * (w_face_n - coeff * (w_face_n - w[i, j, k]))
        end
    end

    # Openings (outlet patches)
    for op in bc_set.openings
        if op.flow_type != OpeningOutlet || !any(isnan, op.velocity)
            continue
        end
        face = op.boundary
        region_check = (i, j, k) -> BoundaryConditions.is_in_opening(op, grid, i, j, k)
        uc_region = BoundaryConditions.compute_region_average_velocity(u, v, w, grid, face, region_check)
        uc_region = clip_uc(uc_region)

        if face == :x_max
            i = mx - 2
            fi = mx - 1
            @inbounds for k in 3:mz-2, j in 3:my-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[fi, j, k]
                    u_face_n = u_face_prev[fi, j, k]
                    buffers.u_face_x[fi, j, k] = m0 * (u_face_n - 2.0 * dt * uc_region * (u_face_n - u[i, j, k]) / dx)
                end
            end
        elseif face == :x_min
            i = 3
            fi = 3
            @inbounds for k in 3:mz-2, j in 3:my-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[fi-1, j, k]
                    u_face_n = u_face_prev[fi, j, k]
                    buffers.u_face_x[fi, j, k] = m0 * (u_face_n - 2.0 * dt * uc_region * (u_face_n - u[i, j, k]) / dx)
                end
            end
        elseif face == :y_max
            j = my - 2
            fj = my - 1
            @inbounds for k in 3:mz-2, i in 3:mx-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[i, fj, k]
                    v_face_n = v_face_prev[i, fj, k]
                    buffers.v_face_y[i, fj, k] = m0 * (v_face_n - 2.0 * dt * uc_region * (v_face_n - v[i, j, k]) / dy)
                end
            end
        elseif face == :y_min
            j = 3
            fj = 3
            @inbounds for k in 3:mz-2, i in 3:mx-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[i, fj-1, k]
                    v_face_n = v_face_prev[i, fj, k]
                    buffers.v_face_y[i, fj, k] = m0 * (v_face_n - 2.0 * dt * uc_region * (v_face_n - v[i, j, k]) / dy)
                end
            end
        elseif face == :z_max
            k = mz - 2
            fk = mz - 1
            @inbounds for j in 3:my-2, i in 3:mx-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[i, j, fk]
                    dz = grid.z_face[k+1] - grid.z_face[k]
                    w_face_n = w_face_prev[i, j, fk]
                    buffers.w_face_z[i, j, fk] = m0 * (w_face_n - 2.0 * dt * uc_region * (w_face_n - w[i, j, k]) / dz)
                end
            end
        elseif face == :z_min
            k = 3
            fk = 3
            @inbounds for j in 3:my-2, i in 3:mx-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[i, j, fk-1]
                    dz = grid.z_face[k+1] - grid.z_face[k]
                    w_face_n = w_face_prev[i, j, fk]
                    buffers.w_face_z[i, j, fk] = m0 * (w_face_n - 2.0 * dt * uc_region * (w_face_n - w[i, j, k]) / dz)
                end
            end
        end
    end
end

"""
    enforce_outflow_face_continuity!(buffers, grid, bc_set)

Outflow境界のセルフェイス速度を連続の式から決定する。
"""
function enforce_outflow_face_continuity!(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    mask = buffers.mask

    if bc_set.x_max.velocity_type == Outflow
        i = mx - 2
        fi = mx - 1
        @inbounds for k in 3:mz-2, j in 3:my-2
            m0 = mask[i, j, k] * mask[fi, j, k]
            dz = grid.z_face[k+1] - grid.z_face[k]
            term_y = (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy
            term_z = (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
            buffers.u_face_x[fi, j, k] = m0 * (buffers.u_face_x[i, j, k] - dx * (term_y + term_z))
        end
    end
    if bc_set.x_min.velocity_type == Outflow
        i = 3
        fi = 3
        @inbounds for k in 3:mz-2, j in 3:my-2
            m0 = mask[i, j, k] * mask[fi-1, j, k]
            dz = grid.z_face[k+1] - grid.z_face[k]
            term_y = (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy
            term_z = (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
            buffers.u_face_x[fi, j, k] = m0 * (buffers.u_face_x[i+1, j, k] + dx * (term_y + term_z))
        end
    end

    if bc_set.y_max.velocity_type == Outflow
        j = my - 2
        fj = my - 1
        @inbounds for k in 3:mz-2, i in 3:mx-2
            m0 = mask[i, j, k] * mask[i, fj, k]
            dz = grid.z_face[k+1] - grid.z_face[k]
            term_x = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx
            term_z = (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
            buffers.v_face_y[i, fj, k] = m0 * (buffers.v_face_y[i, j, k] - dy * (term_x + term_z))
        end
    end
    if bc_set.y_min.velocity_type == Outflow
        j = 3
        fj = 3
        @inbounds for k in 3:mz-2, i in 3:mx-2
            m0 = mask[i, j, k] * mask[i, fj-1, k]
            dz = grid.z_face[k+1] - grid.z_face[k]
            term_x = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx
            term_z = (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
            buffers.v_face_y[i, fj, k] = m0 * (buffers.v_face_y[i, j+1, k] + dy * (term_x + term_z))
        end
    end

    if bc_set.z_max.velocity_type == Outflow
        k = mz - 2
        fk = mz - 1
        @inbounds for j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k] * mask[i, j, fk]
            term_x = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx
            term_y = (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy
            dz = grid.z_face[k+1] - grid.z_face[k]
            buffers.w_face_z[i, j, fk] = m0 * (buffers.w_face_z[i, j, k] - dz * (term_x + term_y))
        end
    end
    if bc_set.z_min.velocity_type == Outflow
        k = 3
        fk = 3
        @inbounds for j in 3:my-2, i in 3:mx-2
            m0 = mask[i, j, k] * mask[i, j, fk-1]
            term_x = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx
            term_y = (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy
            dz = grid.z_face[k+1] - grid.z_face[k]
            buffers.w_face_z[i, j, fk] = m0 * (buffers.w_face_z[i, j, k+1] + dz * (term_x + term_y))
        end
    end

    # Openings (outlet patches)
    for op in bc_set.openings
        if op.flow_type != OpeningOutlet || !any(isnan, op.velocity)
            continue
        end
        face = op.boundary
        region_check = (i, j, k) -> BoundaryConditions.is_in_opening(op, grid, i, j, k)

        if face == :x_max
            i = mx - 2
            fi = mx - 1
            @inbounds for k in 3:mz-2, j in 3:my-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[fi, j, k]
                    dz = grid.z_face[k+1] - grid.z_face[k]
                    term_y = (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy
                    term_z = (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
                    buffers.u_face_x[fi, j, k] = m0 * (buffers.u_face_x[i, j, k] - dx * (term_y + term_z))
                end
            end
        elseif face == :x_min
            i = 3
            fi = 3
            @inbounds for k in 3:mz-2, j in 3:my-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[fi-1, j, k]
                    dz = grid.z_face[k+1] - grid.z_face[k]
                    term_y = (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy
                    term_z = (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
                    buffers.u_face_x[fi, j, k] = m0 * (buffers.u_face_x[i+1, j, k] + dx * (term_y + term_z))
                end
            end
        elseif face == :y_max
            j = my - 2
            fj = my - 1
            @inbounds for k in 3:mz-2, i in 3:mx-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[i, fj, k]
                    dz = grid.z_face[k+1] - grid.z_face[k]
                    term_x = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx
                    term_z = (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
                    buffers.v_face_y[i, fj, k] = m0 * (buffers.v_face_y[i, j, k] - dy * (term_x + term_z))
                end
            end
        elseif face == :y_min
            j = 3
            fj = 3
            @inbounds for k in 3:mz-2, i in 3:mx-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[i, fj-1, k]
                    dz = grid.z_face[k+1] - grid.z_face[k]
                    term_x = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx
                    term_z = (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
                    buffers.v_face_y[i, fj, k] = m0 * (buffers.v_face_y[i, j+1, k] + dy * (term_x + term_z))
                end
            end
        elseif face == :z_max
            k = mz - 2
            fk = mz - 1
            @inbounds for j in 3:my-2, i in 3:mx-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[i, j, fk]
                    term_x = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx
                    term_y = (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy
                    dz = grid.z_face[k+1] - grid.z_face[k]
                    buffers.w_face_z[i, j, fk] = m0 * (buffers.w_face_z[i, j, k] - dz * (term_x + term_y))
                end
            end
        elseif face == :z_min
            k = 3
            fk = 3
            @inbounds for j in 3:my-2, i in 3:mx-2
                if region_check(i, j, k)
                    m0 = mask[i, j, k] * mask[i, j, fk-1]
                    term_x = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx
                    term_y = (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy
                    dz = grid.z_face[k+1] - grid.z_face[k]
                    buffers.w_face_z[i, j, fk] = m0 * (buffers.w_face_z[i, j, k+1] + dz * (term_x + term_y))
                end
            end
        end
    end
end

"""
    compute_divergence!(buffers, grid, dt, mach2, par)

発散計算とポアソン右辺生成。
rhs = (1/dt) * ∇·u* - (M^2/Δt^2) * p^n
"""
function compute_divergence!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    mach2::Float64,
    par::String
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    inv_dt = 1.0 / dt
    alpha = mach2 / (dt * dt)
    
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        m0 = buffers.mask[i, j, k]
        dz = grid.dz[k]
        # Div u* = (u*_{i+1/2} - u*_{i-1/2})/dx + ...
        # u_face_x[i] is left face (i-1/2)?
        # u_face_x[i+1] is right face (i+1/2).
        
        div = (buffers.u_face_x[i+1, j, k] - buffers.u_face_x[i, j, k]) / dx +
              (buffers.v_face_y[i, j+1, k] - buffers.v_face_y[i, j, k]) / dy +
              (buffers.w_face_z[i, j, k+1] - buffers.w_face_z[i, j, k]) / dz
        
        buffers.rhs[i, j, k] = m0 * (div * inv_dt - alpha * buffers.p_prev[i, j, k])
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
        dp_dx = (p[i, j, k] - p[i-1, j, k]) / dx * mask[i-1, j, k] * mask[i, j, k]
        buffers.u_face_x[i, j, k] -= dt * dp_dx
    end
    
    # Y方向フェイス (j-1/2位置)
    @inbounds for k in 3:mz-2, j in 3:my-1, i in 3:mx-2
        dp_dy = (p[i, j, k] - p[i, j-1, k]) / dy * mask[i, j-1, k] * mask[i, j, k]
        buffers.v_face_y[i, j, k] -= dt * dp_dy
    end
    
    # Z方向フェイス (k-1/2位置)
    @inbounds for k in 3:mz-1, j in 3:my-2, i in 3:mx-2
        dz = grid.z_face[k] - grid.z_face[k-1]
        dp_dz = (p[i, j, k] - p[i, j, k-1]) / dz * mask[i, j, k-1] * mask[i, j, k]
        buffers.w_face_z[i, j, k] -= dt * dp_dz
    end
    
    # 2. セルセンター速度の修正（両側フェイスの圧力勾配の平均）
    @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
        # X方向: 左右フェイスの圧力勾配を平均
        m0 = mask[i, j, k]

        # 左フェイス (i-1/2) の圧力勾配（マスクでNeumann条件）
        dp_dx_left = m0 * mask[i-1, j, k] * (p[i, j, k] - p[i-1, j, k]) / dx
        # 右フェイス (i+1/2) の圧力勾配
        dp_dx_right = m0 * mask[i+1, j, k] * (p[i+1, j, k] - p[i, j, k]) / dx
        
        # 平均勾配で修正
        dp_dx_avg = 0.5 * (dp_dx_left + dp_dx_right)
        buffers.u[i, j, k] = buffers.u_star[i, j, k] - dt * dp_dx_avg
        
        # Y方向
        dp_dy_front = m0 * mask[i, j-1, k] * (p[i, j, k] - p[i, j-1, k]) / dy
        dp_dy_back = m0 * mask[i, j+1, k] * (p[i, j+1, k] - p[i, j, k]) / dy
        
        dp_dy_avg = 0.5 * (dp_dy_front + dp_dy_back)
        buffers.v[i, j, k] = buffers.v_star[i, j, k] - dt * dp_dy_avg
        
        # Z方向
        dz_bottom = grid.z_face[k] - grid.z_face[k-1]
        dz_top = grid.z_face[k+1] - grid.z_face[k]
        
        dp_dz_bottom = m0 * mask[i, j, k-1] * (p[i, j, k] - p[i, j, k-1]) / dz_bottom
        dp_dz_top = m0 * mask[i, j, k+1] * (p[i, j, k+1] - p[i, j, k]) / dz_top
        
        dp_dz_avg = 0.5 * (dp_dz_bottom + dp_dz_top)
        buffers.w[i, j, k] = buffers.w_star[i, j, k] - dt * dp_dz_avg
    end
end

"""
    fractional_step!(buffers, grid, dt, bc_set, poisson_config, Cs, nu, reverse_flow_stabilization, par)

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
    reverse_flow_stabilization::Bool,
    par::String
)::Tuple{Int, Float64}

    # 1. Update Turbulence (nu_t)
    compute_turbulent_viscosity!(buffers, grid, Cs, nu, par)
    
    # 2. Pseudo Velocity
    compute_pseudo_velocity!(buffers, grid, dt, bc_set, par)
    
    # 2.5 Apply BCs to Pseudo Velocity
    apply_velocity_cc_bcs!(
        buffers.u_star, buffers.v_star, buffers.w_star,
        grid, buffers.mask, bc_set, dt;
        u_ref=buffers.u, v_ref=buffers.v, w_ref=buffers.w,
        reverse_flow_stabilization=reverse_flow_stabilization
    )

    # 3. Interpolate to Faces
    copyto!(buffers.flux_u, buffers.u_face_x)
    copyto!(buffers.flux_v, buffers.v_face_y)
    copyto!(buffers.flux_w, buffers.w_face_z)
    interpolate_to_faces!(buffers, grid, par)

    apply_outflow_face_star!(
        buffers, grid, bc_set, dt,
        buffers.flux_u, buffers.flux_v, buffers.flux_w,
        reverse_flow_stabilization
    )
    
    # 4. Divergence / RHS , copy from p to p_prev
    copyto!(buffers.p_prev, buffers.p)
    compute_divergence!(buffers, grid, dt, poisson_config.mach2, par)
    
    # 5. Poisson Solve
    # Temporarily set Inflow/Outflow/Opening mask to 0.0 to enforce Neumann pressure BC (dp/dn=0)
    update_boundary_mask!(buffers.mask, grid, bc_set, 0.0)
    
    alpha = poisson_config.mach2 / (dt * dt)
    converged, iter, res = solve_poisson!(buffers, grid, poisson_config, bc_set, par, alpha)
    
    # Restore Inflow/Outflow/Opening mask to 1.0 for Convection/Velocity update
    update_boundary_mask!(buffers.mask, grid, bc_set, 1.0)
    
    # 6. Correct Velocity (Face and Center)
    correct_velocity!(buffers, grid, dt, par)
    
    # 7. Apply BCs to Face Velocities
    # セルフェイス速度への境界条件適用（周期境界、対称境界など）
    apply_velocity_cf_bcs!(buffers.u_face_x, buffers.v_face_y, buffers.w_face_z, grid, buffers.mask, bc_set, dt)
    enforce_outflow_face_continuity!(buffers, grid, bc_set)
    
    # 8. Apply BCs to Center Velocities
    apply_velocity_cc_bcs!(
        buffers.u, buffers.v, buffers.w,
        grid, buffers.mask, bc_set, dt;
        reverse_flow_stabilization=reverse_flow_stabilization
    )
    
    return (iter, res)
end

end # module FractionalStep
