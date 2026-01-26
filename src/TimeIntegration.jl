module TimeIntegration

using ..Common
using ..Fields
using ..Grid
using ..BoundaryConditions
using ..PressureSolver
using ..Turbulence
using ..Convection
using ..Diffusion
using ..FractionalStep: interpolate_to_faces!, compute_divergence!, correct_velocity!, compute_pseudo_velocity!, fractional_step!

export TimeScheme, TimeConfig, advance!, compute_dt
export Euler, RK2, RK4

@enum TimeScheme begin
    Euler
    RK2
    RK4
end

struct TimeConfig
    scheme::TimeScheme
    Co::Float64
    dt_fixed::Float64
end

"""
    compute_dt(buffers, grid, Co, nu)

時間刻みを計算（拡散数チェックは別処理）。
"""
function compute_dt(
    buffers::CFDBuffers,
    grid::GridData,
    Co::Float64,
    nu::Float64
)::Float64
    u_max = 0.0
    v_max = 0.0
    w_max = 0.0
    
    # Parallel reduction capable logic
    @inbounds for i in eachindex(buffers.u)
        u_max = max(u_max, abs(buffers.u[i]))
        v_max = max(v_max, abs(buffers.v[i]))
        w_max = max(w_max, abs(buffers.w[i]))
    end
    
    U_max = sqrt(u_max^2 + v_max^2 + w_max^2)
    # Use max(U_max, 1.0) to ensure stable time step at initialization
    # U*_ref = max(U*_max, 1.0) where 1.0 is the non-dimensional reference velocity
    U_ref = max(U_max, 1.0)
    
    dx_min = min(grid.dx, grid.dy)
    if !isempty(grid.dz)
        dx_min = min(dx_min, minimum(grid.dz))
    end
    
    dt_adv = Co * dx_min / U_ref
    
    return dt_adv
end

"""
    compute_dt(grid, Co, nu, U_ref)

時間刻みを計算（初期速度などからU_refを指定）。
"""
function compute_dt(
    grid::GridData,
    Co::Float64,
    nu::Float64,
    U_ref::Float64
)::Float64
    U_ref = max(U_ref, 1.0)
    dx_min = min(grid.dx, grid.dy)
    if !isempty(grid.dz)
        dx_min = min(dx_min, minimum(grid.dz))
    end

    return Co * dx_min / U_ref
end

# --- Helper: Compute Flux H(u) ---
function compute_flux!(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    Cs::Float64,
    nu::Float64,
    par::String
)
    # 1. Update Turbulence (nu_t depends on current Velocity u)
    compute_turbulent_viscosity!(buffers.nu_t, buffers, grid, Cs, nu, par)
    
    # 2. Clear Flux
    fill!(buffers.flux_u, 0.0)
    fill!(buffers.flux_v, 0.0)
    fill!(buffers.flux_w, 0.0)
    
    # 3. Add Convection & Diffusion
    add_convection_flux!(buffers, grid, bc_set, par)
    add_diffusion_flux!(buffers, grid, bc_set, par)
    
    # buffers.flux_u/v/w now holds H(u)
end

# --- Helper: Projection Step (u* -> u) ---
function projection_step!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    bc_set::BoundaryConditionSet,
    poisson_config::PoissonConfig,
    par::String;
    u_ref::Array{Float64, 3}=buffers.u,
    v_ref::Array{Float64, 3}=buffers.v,
    w_ref::Array{Float64, 3}=buffers.w
)::Tuple{Int, Float64}
    # 0. Apply BCs to u*
    apply_velocity_bcs!(
        buffers.u_star, buffers.v_star, buffers.w_star,
        grid, buffers.mask, bc_set, dt;
        u_ref=u_ref, v_ref=v_ref, w_ref=w_ref
    )

    # 1. Interpolate to faces
    interpolate_to_faces!(buffers, grid, par)
    
    # 2. Divergence
    copyto!(buffers.p_prev, buffers.p)
    compute_divergence!(buffers, grid, dt, poisson_config.mach2, par)
    
    # 3. Poisson
    alpha = poisson_config.mach2 / (dt * dt)
    converged, iter, res = solve_poisson!(buffers, grid, poisson_config, bc_set, par, alpha)
    
    # 4. Correct
    correct_velocity!(buffers, grid, dt, par)
    
    # 5. Apply BCs
    apply_boundary_conditions!(buffers, grid, bc_set, dt, par)
    
    return (iter, res)
end


function euler_step!(
    buffers::CFDBuffers,
    grid::GridData,
    dt::Float64,
    bc_set::BoundaryConditionSet,
    poisson_config::PoissonConfig,
    Cs::Float64,
    nu::Float64,
    par::String
)
    # Reuse Full Fractional Step
    return fractional_step!(buffers, grid, dt, bc_set, poisson_config, Cs, nu, par)
end

function rk2_step!(
    buffers::CFDBuffers,
    rk_buffers::RKBuffers,
    grid::GridData,
    dt::Float64,
    bc_set::BoundaryConditionSet,
    poisson_config::PoissonConfig,
    Cs::Float64,
    nu::Float64,
    par::String
)
    # Midpoint Method (Predictor-Corrector)
    # 1. Store u^n -> rk1
    copy!(rk_buffers.u_rk1, buffers.u)
    copy!(rk_buffers.v_rk1, buffers.v)
    copy!(rk_buffers.w_rk1, buffers.w)
    
    # 2. Stage 1: Predict u^{n+1/2} (Euler with dt/2)
    # We can use fractional_step!(dt/2) directly on buffers.u
    fractional_step!(buffers, grid, dt * 0.5, bc_set, poisson_config, Cs, nu, par)
    # buffers.u is now u^{n+1/2} (approx)
    
    # 3. Stage 2: H(u^{n+1/2})
    # stored in buffers.flux_u/v/w
    compute_flux!(buffers, grid, bc_set, Cs, nu, par)
    # flux is H(u^{n+1/2})
    
    # 4. Update u* = u^n + dt * flux
    # Restore u^n from rk1 for base
    # But add flux to it
    @inbounds @. buffers.u_star = rk_buffers.u_rk1 + dt * buffers.flux_u
    @inbounds @. buffers.v_star = rk_buffers.v_rk1 + dt * buffers.flux_v
    @inbounds @. buffers.w_star = rk_buffers.w_rk1 + dt * buffers.flux_w
    
    # 5. Projection u* -> u^{n+1}
    iter, res = projection_step!(
        buffers, grid, dt, bc_set, poisson_config, par;
        u_ref=rk_buffers.u_rk1, v_ref=rk_buffers.v_rk1, w_ref=rk_buffers.w_rk1
    )
    
    return (iter, res)
end

function rk4_step!(
    buffers::CFDBuffers,
    rk_buffers::RKBuffers,
    grid::GridData,
    dt::Float64,
    bc_set::BoundaryConditionSet,
    poisson_config::PoissonConfig,
    Cs::Float64,
    nu::Float64,
    par::String
)
    # 4-Stage RK with Projection
    # Save u^n
    copy!(rk_buffers.u_rk1, buffers.u) # reusing rk1 as storage for u^n ??? NO.
    # rk1..rk4 logic:
    # We store k1..k4 in rk1..rk4?
    # We need u^n preserved. 
    # Let's say u^n is stored in `u_rk4` temporarily? No, rk4 is needed for k4.
    # We need a dedicated buffer for u^n.
    # Available: u_rk1, u_rk2, u_rk3, u_rk4.
    # Total 4 slots.
    # Algorithm needs: u^n, k1, k2, k3, k4. (5 vectors).
    # But we can accumulate result: u_{acc} = u^n + ...
    # And we need u^n for each stage base.
    # So we MUST keep u^n.
    # Use `u_rk1` for u^n.
    # Use `u_rk2` for k1, `u_rk3` for k2, `u_rk4` for k3.
    # What about k4? We can compute and add immediately.
    # But final sum needs k1, k2, k3.
    # 4 slots is enough if careful.
    
    # Alias
    u_n = (rk_buffers.u_rk1, rk_buffers.v_rk1, rk_buffers.w_rk1)
    k1  = (rk_buffers.u_rk2, rk_buffers.v_rk2, rk_buffers.w_rk2)
    k2  = (rk_buffers.u_rk3, rk_buffers.v_rk3, rk_buffers.w_rk3)
    k3  = (rk_buffers.u_rk4, rk_buffers.v_rk4, rk_buffers.w_rk4)
    # k4 uses flux buffers directly
    
    # 1. Save u^n
    copy!(u_n[1], buffers.u)
    copy!(u_n[2], buffers.v)
    copy!(u_n[3], buffers.w)
    
    # --- Stage 1 ---
    # Calc k1 = H(u^n)
    # stored in buffers.flux_u/v/w
    compute_flux!(buffers, grid, bc_set, Cs, nu, par)
    copy!(k1[1], buffers.flux_u)
    copy!(k1[2], buffers.flux_v)
    copy!(k1[3], buffers.flux_w)
    
    # Update u for Stage 2: u1 = u^n + 0.5 * dt * k1
    @inbounds @. buffers.u_star = u_n[1] + 0.5 * dt * k1[1]
    @inbounds @. buffers.v_star = u_n[2] + 0.5 * dt * k1[2]
    @inbounds @. buffers.w_star = u_n[3] + 0.5 * dt * k1[3]
    
    # Project u1 (inplace to buffers.u)
    projection_step!(
        buffers, grid, 0.5*dt, bc_set, poisson_config, par;
        u_ref=u_n[1], v_ref=u_n[2], w_ref=u_n[3]
    )
    
    # --- Stage 2 ---
    # Calc k2 = H(u1)
    # stored in buffers.flux_u/v/w
    compute_flux!(buffers, grid, bc_set, Cs, nu, par)
    copy!(k2[1], buffers.flux_u)
    copy!(k2[2], buffers.flux_v)
    copy!(k2[3], buffers.flux_w)
    
    # Update u for Stage 3: u2 = u^n + 0.5 * dt * k2
    @inbounds @. buffers.u_star = u_n[1] + 0.5 * dt * k2[1]
    @inbounds @. buffers.v_star = u_n[2] + 0.5 * dt * k2[2]
    @inbounds @. buffers.w_star = u_n[3] + 0.5 * dt * k2[3]
    
    projection_step!(
        buffers, grid, 0.5*dt, bc_set, poisson_config, par;
        u_ref=u_n[1], v_ref=u_n[2], w_ref=u_n[3]
    )
    
    # --- Stage 3 ---
    # Calc k3 = H(u2)
    # stored in buffers.flux_u/v/w
    compute_flux!(buffers, grid, bc_set, Cs, nu, par)
    copy!(k3[1], buffers.flux_u)
    copy!(k3[2], buffers.flux_v)
    copy!(k3[3], buffers.flux_w)
    
    # Update u for Stage 4: u3 = u^n + dt * k3
    @inbounds @. buffers.u_star = u_n[1] + dt * k3[1]
    @inbounds @. buffers.v_star = u_n[2] + dt * k3[2]
    @inbounds @. buffers.w_star = u_n[3] + dt * k3[3]
    
    projection_step!(
        buffers, grid, dt, bc_set, poisson_config, par;
        u_ref=u_n[1], v_ref=u_n[2], w_ref=u_n[3]
    )
    
    # Calc k4 = H(u3)
    # stored in buffers.flux_u/v/w
    compute_flux!(buffers, grid, bc_set, Cs, nu, par)
    
    # Final Accumulation
    # u* = u^n + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    # using buffers.u_star as target
    dt6 = dt / 6.0
    @inbounds @. buffers.u_star = u_n[1] + dt6 * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + buffers.flux_u)
    @inbounds @. buffers.v_star = u_n[2] + dt6 * (k1[2] + 2.0*k2[2] + 2.0*k3[2] + buffers.flux_v)
    @inbounds @. buffers.w_star = u_n[3] + dt6 * (k1[3] + 2.0*k2[3] + 2.0*k3[3] + buffers.flux_w)
    
    # Final Projection
    iter, res = projection_step!(
        buffers, grid, dt, bc_set, poisson_config, par;
        u_ref=u_n[1], v_ref=u_n[2], w_ref=u_n[3]
    )
    
    return (iter, res)
end

"""
    advance!(buffers, grid, bc_set, time_config, poisson_config, Cs, par; rk_buffers=nothing)

1タイムステップ進行。
"""
# === Main Entry Point: Advance One Time Step ===
function advance!(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    time_config::TimeConfig,
    poisson_config::PoissonConfig,
    Cs::Float64,
    nu::Float64,
    dt_fixed::Float64,  # Fixed time step (calculated once at initialization)
    par::String;
    rk_buffers::Union{RKBuffers, Nothing} = nothing
)::Tuple{Int, Float64}
    # Returns: (圧力反復回数, 圧力残差)
    
    if time_config.scheme == Euler
        return euler_step!(buffers, grid, dt_fixed, bc_set, poisson_config, Cs, nu, par)
    elseif time_config.scheme == RK2
        if isnothing(rk_buffers)
            error("RK2 selected but rk_buffers not provided.")
        end
        return rk2_step!(buffers, rk_buffers, grid, dt_fixed, bc_set, poisson_config, Cs, nu, par)
    elseif time_config.scheme == RK4
        if isnothing(rk_buffers)
            error("RK4 selected but rk_buffers not provided.")
        end
        return rk4_step!(buffers, rk_buffers, grid, dt_fixed, bc_set, poisson_config, Cs, nu, par)
    else
        error("Unknown time integration scheme: $(time_config.scheme)")
    end
end

end # module TimeIntegration
