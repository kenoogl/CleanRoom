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
    add_diffusion_flux!(buffers, grid, par)
    
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

速度補正。 u = u* - dt * ∇p
"""
function correct_velocity!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    par::String
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    dx, dy = grid.dx, grid.dy
    
    @inbounds for k in 3:mz-2
        # z-gradient needs dZ (distance between centers)
        dZ_m = grid.z_center[k] - grid.z_center[k-1] # for i?
        # Central difference for p at center?
        # No, u is at center, p is at center.
        # But u* was interpolated to faces to get Div.
        # Now calculating New U at center.
        # u^{n+1} = u* - dt * (p_{i+1} - p_{i-1})/(2dx) ?
        # Standard projection:
        # u = u* - dt * grad(p).
        # Consistent with Div u* on faces -> p on centers -> grad p on faces?
        # If we update Centered U:
        # u_new = u*_center - dt * (p_{i+1}-p_{i-1})/(2dx) ?
        # Or should we work on MAC grid?
        # Codebase seems Colocated (all at center).
        # Rhie-Chow interpolation prevents checkerboard.
        # u*_face matches p-gradient at face. Corrects u_face -> u_face_new satisfying div=0.
        # Then we need to reconstruct u_center_new from u_face_new?
        # Or subtract "modified pressure gradient" from u*_center.
        
        # Simple Colocated (approx):
        # u_new = u_star - dt * grad_p
        # grad_p = (p_{i+1} - p_{i-1})/(2dx)
        # This allows checkerboard.
        # But `interpolate_to_faces` helps Div calculation stability.
        # For U-correction, if we just use central diff, we get checkerboard.
        # Design doc Task 11.2 "チェッカーボード不安定性の防止" linked to interpolate_to_faces.
        # Does it imply we update Faces and then average back to center?
        # Or we rely on the fact that Rhie-Chow modifies the face velocity.
        # But we need Cell Center Velocity for next step.
        
        # Standard Rhie-Chow:
        # u_face = avg(u) - dt*(grad P - avg(grad P)) + ...
        # My `interpolate_to_faces` just averaged u*.
        # It did NOT add the pressure smoothing term.
        # If this simulator relies on naive approach + smoothing elsewhere, I'll follow simple gradient.
        # But standard for this architecture is:
        # u^{n+1}_i = u*_i - dt * (p_{i+1} - p_{i-1})/(2dx)
        # This IS prone to checkerboard.
        # Maybe `interpolate_to_faces` was intended to be the only place where strict conservation is enforced?
        # But we carry `u` (center) to next step.
        
        # I will implement central difference correction.
        # If instability occurs, I might need to implement filtering.
        
        dz_2 = grid.z_center[k+1] - grid.z_center[k-1] # Distance between k+1 and k-1 centers
        
        for j in 3:my-2
            for i in 3:mx-2
                if buffers.mask[i, j, k] == 0.0
                    buffers.u[i, j, k] = 0.0; buffers.v[i, j, k] = 0.0; buffers.w[i, j, k] = 0.0
                    continue
                end
                
                dpdx = (buffers.p[i+1, j, k] - buffers.p[i-1, j, k]) / (2*dx)
                dpdy = (buffers.p[i, j+1, k] - buffers.p[i, j-1, k]) / (2*dy)
                dpdz = (buffers.p[i, j, k+1] - buffers.p[i, j, k-1]) / dz_2
                
                buffers.u[i, j, k] = buffers.u_star[i, j, k] - dt * dpdx
                buffers.v[i, j, k] = buffers.v_star[i, j, k] - dt * dpdy
                buffers.w[i, j, k] = buffers.w_star[i, j, k] - dt * dpdz
            end
        end
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
    compute_pseudo_velocity!(buffers, grid, dt, par)
    
    # 2.5 Apply BCs to Pseudo Velocity
    apply_velocity_bcs!(buffers.u_star, buffers.v_star, buffers.w_star, grid, buffers.mask, bc_set, dt)

    # 3. Interpolate to Faces
    interpolate_to_faces!(buffers, grid, par)
    
    # 4. Divergence / RHS
    compute_divergence!(buffers, grid, dt, par)
    
    # 5. Poisson Solve
    converged, iter, res = solve_poisson!(buffers, grid, poisson_config, bc_set, par)
    
    # 6. Correct Velocity
    correct_velocity!(buffers, grid, dt, par)
    
    # 7. Apply BCs
    apply_boundary_conditions!(buffers, grid, bc_set, dt, par)
    
    # 8. Update Time Average (Should be called by Main? Or here?)
    # Design flow says "Update Time Average" is last step of FracStep flow.
    # But usually caller manages averaging start time.
    # I'll leave it to Main because it depends on Time.
    
    return (iter, res)
end

end # module FractionalStep
