module BoundaryConditions

using ..Common
using ..Grid
using ..Fields

export VelocityBCType, ExternalBC, InletOutlet, InternalBoundary, BoundaryConditionSet
export Dirichlet, Neumann, Outflow, Periodic
export apply_boundary_conditions!, apply_outflow!, apply_velocity_bcs!
export apply_periodic_velocity!, apply_periodic_pressure!

@enum VelocityBCType begin
    Dirichlet
    Neumann
    Outflow
    Periodic
end

struct ExternalBC
    velocity_type::VelocityBCType
    velocity_value::NTuple{3, Float64}
end

struct InletOutlet
    type::Symbol
    position::NTuple{3, Float64}
    size::NTuple{2, Float64}
    normal::NTuple{3, Int}
    condition::VelocityBCType
    velocity::NTuple{3, Float64}
end

struct InternalBoundary
    type::Symbol
    region_min::NTuple{3, Float64}
    region_max::NTuple{3, Float64}
    center::NTuple{3, Float64}
    radius::Float64
    height::Float64
    axis::Symbol
    normal::NTuple{3, Int}
    velocity::NTuple{3, Float64}
end

struct BoundaryConditionSet
    x_min::ExternalBC
    x_max::ExternalBC
    y_min::ExternalBC
    y_max::ExternalBC
    z_min::ExternalBC
    z_max::ExternalBC
    inlets::Vector{InletOutlet}
    outlets::Vector{InletOutlet}
    internal_boundaries::Vector{InternalBoundary}
end

"""
    compute_face_average_velocity(buffers, grid, face)

境界面の平均流速を計算する。
"""
function compute_face_average_velocity(
    u::Array{Float64, 3},
    v::Array{Float64, 3},
    w::Array{Float64, 3},
    grid::GridData,
    face::Symbol
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    sum_vel = 0.0
    count = 0
    
    if face == :x_min
        i = 3 # Boundary cell (Real start)
        for k in 3:mz-2, j in 3:my-2
            sum_vel += u[i, j, k]
            count += 1
        end
    elseif face == :x_max
        i = mx-2
        for k in 3:mz-2, j in 3:my-2
            sum_vel += u[i, j, k]
            count += 1
        end
    elseif face == :y_min
        j = 3
        for k in 3:mz-2, i in 3:mx-2
            sum_vel += v[i, j, k]
            count += 1
        end
    elseif face == :y_max
        j = my-2
        for k in 3:mz-2, i in 3:mx-2
            sum_vel += v[i, j, k]
            count += 1
        end
    elseif face == :z_min
        k = 3
        for j in 3:my-2, i in 3:mx-2
            sum_vel += w[i, j, k]
            count += 1
        end
    elseif face == :z_max
        k = mz-2
        for j in 3:my-2, i in 3:mx-2
            sum_vel += w[i, j, k]
            count += 1
        end
    end
    
    return (count > 0) ? sum_vel / count : 0.0
end

"""
    apply_outflow!(phi, grid, face, dt, Uc)

対流流出条件を適用する。
"""
function apply_outflow!(
    phi::Array{Float64, 3},
    grid::GridData,
    face::Symbol,
    dt::Float64,
    Uc::Float64
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    # 1st order upwind: phi_new = phi - Uc * dt/dx * (phi - phi_inner)
    # If Uc > 0 (outflow), we use values from upstream/inner.
    # Boundary is at ghost cell? No, usually boundary condition sets values at Ghost Cells
    # OR at boundary face.
    # Here Finite Volume with Ghost Cells.
    # Dirichlet: Set Ghost = 2*Val - Inner (if Val is at face) or Ghost=Val (if Val is at center of ghost).
    # Outflow: Set Ghost such that dphi/dt + U dphi/dn = 0.
    
    # x_max (i=mx-1, mx-0 ghosts). Inner is i=mx-2.
    # Outflow at interface i+1/2 (between mx-2 and mx-1).
    # We update value at Ghost Cell i=mx-1.
    # dphi/dt + U (phi_{mx-1} - phi_{mx-2})/dx = 0?
    # phi_{mx-1}^{n+1} = phi_{mx-1}^n - U*dt/dx * (phi_{mx-1}^n - phi_{mx-2}^n).
    
    # Need to verify if we update Ghost cells "evolutively" or just extrapolation from Inner.
    # Simple extrapolation (zero gradient) is Neumann.
    # Convective outflow solves equation.
    
    if face == :x_max
        # i=mx-1 is ghost 1.
        dx = grid.dx
        c = Uc * dt / dx
        @inbounds for k in 1:mz, j in 1:my
            phi[mx-1, j, k] = phi[mx-1, j, k] - c * (phi[mx-1, j, k] - phi[mx-2, j, k])
            phi[mx, j, k] = phi[mx-1, j, k] # Copy to outer ghost
        end
    elseif face == :x_min
        # i=2 is ghost. Inner i=3. Flow is -x direction (Uc < 0).
        # Eq: dphi/dt + U dphi/dx = 0.
        # U is negative. dphi/dx approx (phi_3 - phi_2)/dx.
        # phi_2_new = phi_2 - U*dt * (phi_3 - phi_2)/dx
        # Since U approx -abs(U), c = U*dt/dx is negative.
        # Term -c * (...) is positive contribution from inside.
        dx = grid.dx
        c = Uc * dt / dx
        @inbounds for k in 1:mz, j in 1:my
            phi[2, j, k] = phi[2, j, k] - c * (phi[3, j, k] - phi[2, j, k])
            phi[1, j, k] = phi[2, j, k]
        end
    elseif face == :y_max
        dy = grid.dy
        c = Uc * dt / dy
        @inbounds for k in 1:mz, i in 1:mx
            phi[i, my-1, k] = phi[i, my-1, k] - c * (phi[i, my-1, k] - phi[i, my-2, k])
            phi[i, my, k] = phi[i, my-1, k]
        end
    elseif face == :y_min
        dy = grid.dy
        c = Uc * dt / dy
        @inbounds for k in 1:mz, i in 1:mx
            phi[i, 2, k] = phi[i, 2, k] - c * (phi[i, 3, k] - phi[i, 2, k])
            phi[i, 1, k] = phi[i, 2, k]
        end
    elseif face == :z_max
        # dz at boundary?
        # Ghost cell width? Assume grid.dz[mz-1] ?
        # Or Just use inner cell width.
        dz = grid.dz[mz-2]
        c = Uc * dt / dz
        @inbounds for j in 1:my, i in 1:mx
            phi[i, j, mz-1] = phi[i, j, mz-1] - c * (phi[i, j, mz-1] - phi[i, j, mz-2])
            phi[i, j, mz] = phi[i, j, mz-1]
        end
    elseif face == :z_min
        dz = grid.dz[3]
        c = Uc * dt / dz
        @inbounds for j in 1:my, i in 1:mx
            phi[i, j, 2] = phi[i, j, 2] - c * (phi[i, j, 3] - phi[i, j, 2])
            phi[i, j, 1] = phi[i, j, 2]
        end
    end
end

"""
    set_boundary_value!(phi, grid, face, val, type)

Helper to set Dirichlet/Neumann on a face.
"""
function set_boundary_value!(
    phi::Array{Float64, 3},
    grid::GridData,
    face::Symbol,
    val::Float64,
    type::VelocityBCType
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    if type == Dirichlet
        if face == :x_min
             @inbounds phi[1:2, :, :] .= val
        elseif face == :x_max
             @inbounds phi[mx-1:mx, :, :] .= val
        elseif face == :y_min
             @inbounds phi[:, 1:2, :] .= val
        elseif face == :y_max
             @inbounds phi[:, my-1:my, :] .= val
        elseif face == :z_min
             @inbounds phi[:, :, 1:2] .= val
        elseif face == :z_max
             @inbounds phi[:, :, mz-1:mz] .= val
        end
    elseif type == Neumann
        if face == :x_min
             @inbounds @. phi[2, :, :] = phi[3, :, :]; phi[1, :, :] = phi[2, :, :]
        elseif face == :x_max
             @inbounds @. phi[mx-1, :, :] = phi[mx-2, :, :]; phi[mx, :, :] = phi[mx-1, :, :]
        elseif face == :y_min
             @inbounds @. phi[:, 2, :] = phi[:, 3, :]; phi[:, 1, :] = phi[:, 2, :]
        elseif face == :y_max
             @inbounds @. phi[:, my-1, :] = phi[:, my-2, :]; phi[:, my, :] = phi[:, my-1, :]
        elseif face == :z_min
             @inbounds @. phi[:, :, 2] = phi[:, :, 3]; phi[:, :, 1] = phi[:, :, 2]
        elseif face == :z_max
             @inbounds @. phi[:, :, mz-1] = phi[:, :, mz-2]; phi[:, :, mz] = phi[:, :, mz-1]
        end
    end
end


"""
    apply_periodic_velocity!(u, v, w, grid, axis)

速度成分の周期境界条件を適用する。
"""
function apply_periodic_velocity!(
    u::Array{Float64, 3},
    v::Array{Float64, 3},
    w::Array{Float64, 3},
    grid::GridData,
    axis::Symbol
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    if axis == :y
        @inbounds for k in 1:mz, i in 1:mx
            u[i, 1, k] = u[i, my-3, k]
            u[i, 2, k] = u[i, my-2, k]
            u[i, my-1, k] = u[i, 3, k]
            u[i, my, k] = u[i, 4, k]
            
            v[i, 1, k] = v[i, my-3, k]
            v[i, 2, k] = v[i, my-2, k]
            v[i, my-1, k] = v[i, 3, k]
            v[i, my, k] = v[i, 4, k]
            
            w[i, 1, k] = w[i, my-3, k]
            w[i, 2, k] = w[i, my-2, k]
            w[i, my-1, k] = w[i, 3, k]
            w[i, my, k] = w[i, 4, k]
        end
    elseif axis == :x
        @inbounds for k in 1:mz, j in 1:my
            u[1, j, k] = u[mx-3, j, k]
            u[2, j, k] = u[mx-2, j, k]
            u[mx-1, j, k] = u[3, j, k]
            u[mx, j, k] = u[4, j, k]
            
            v[1, j, k] = v[mx-3, j, k]
            v[2, j, k] = v[mx-2, j, k]
            v[mx-1, j, k] = v[3, j, k]
            v[mx, j, k] = v[4, j, k]
            
            w[1, j, k] = w[mx-3, j, k]
            w[2, j, k] = w[mx-2, j, k]
            w[mx-1, j, k] = w[3, j, k]
            w[mx, j, k] = w[4, j, k]
        end
    elseif axis == :z
        @inbounds for j in 1:my, i in 1:mx
            u[i, j, 1] = u[i, j, mz-3]
            u[i, j, 2] = u[i, j, mz-2]
            u[i, j, mz-1] = u[i, j, 3]
            u[i, j, mz] = u[i, j, 4]
            
            v[i, j, 1] = v[i, j, mz-3]
            v[i, j, 2] = v[i, j, mz-2]
            v[i, j, mz-1] = v[i, j, 3]
            v[i, j, mz] = v[i, j, 4]
            
            w[i, j, 1] = w[i, j, mz-3]
            w[i, j, 2] = w[i, j, mz-2]
            w[i, j, mz-1] = w[i, j, 3]
            w[i, j, mz] = w[i, j, 4]
        end
    end
end

"""
    apply_periodic_pressure!(p, grid, axis)

圧力の周期境界条件を適用する。
"""
function apply_periodic_pressure!(
    p::Array{Float64, 3},
    grid::GridData,
    axis::Symbol
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    if axis == :y
        @inbounds for k in 1:mz, i in 1:mx
            p[i, 1, k] = p[i, my-3, k]
            p[i, 2, k] = p[i, my-2, k]
            p[i, my-1, k] = p[i, 3, k]
            p[i, my, k] = p[i, 4, k]
        end
    elseif axis == :x
        @inbounds for k in 1:mz, j in 1:my
            p[1, j, k] = p[mx-3, j, k]
            p[2, j, k] = p[mx-2, j, k]
            p[mx-1, j, k] = p[3, j, k]
            p[mx, j, k] = p[4, j, k]
        end
    elseif axis == :z
        @inbounds for j in 1:my, i in 1:mx
            p[i, j, 1] = p[i, j, mz-3]
            p[i, j, 2] = p[i, j, mz-2]
            p[i, j, mz-1] = p[i, j, 3]
            p[i, j, mz] = p[i, j, 4]
        end
    end
end


"""
    apply_velocity_bcs!(u, v, w, grid, bc_set, dt)

指定された速度場(u, v, w)に対して境界条件を適用する。
"""
function apply_velocity_bcs!(
    u::Array{Float64, 3},
    v::Array{Float64, 3},
    w::Array{Float64, 3},
    grid::GridData,
    bc_set::BoundaryConditionSet,
    dt::Float64
)
    # 1. External Boundaries
    function apply_ext(bc, face)
        if bc.velocity_type == Dirichlet
             set_boundary_value!(u, grid, face, bc.velocity_value[1], Dirichlet)
             set_boundary_value!(v, grid, face, bc.velocity_value[2], Dirichlet)
             set_boundary_value!(w, grid, face, bc.velocity_value[3], Dirichlet)
        elseif bc.velocity_type == Neumann
             set_boundary_value!(u, grid, face, 0.0, Neumann)
             set_boundary_value!(v, grid, face, 0.0, Neumann)
             set_boundary_value!(w, grid, face, 0.0, Neumann)
        elseif bc.velocity_type == Outflow
             Uc = compute_face_average_velocity(u, v, w, grid, face)
             apply_outflow!(u, grid, face, dt, Uc)
             apply_outflow!(v, grid, face, dt, Uc)
             apply_outflow!(w, grid, face, dt, Uc)
        end
    end

    # Check for periodic BC pairs
    y_periodic = bc_set.y_min.velocity_type == Periodic && bc_set.y_max.velocity_type == Periodic
    x_periodic = bc_set.x_min.velocity_type == Periodic && bc_set.x_max.velocity_type == Periodic
    z_periodic = bc_set.z_min.velocity_type == Periodic && bc_set.z_max.velocity_type == Periodic
    
    if y_periodic
        apply_periodic_velocity!(u, v, w, grid, :y)
    else
        apply_ext(bc_set.y_min, :y_min)
        apply_ext(bc_set.y_max, :y_max)
    end
    
    if x_periodic
        apply_periodic_velocity!(u, v, w, grid, :x)
    else
        apply_ext(bc_set.x_min, :x_min)
        apply_ext(bc_set.x_max, :x_max)
    end
    
    if z_periodic
        apply_periodic_velocity!(u, v, w, grid, :z)
    else
        apply_ext(bc_set.z_min, :z_min)
        apply_ext(bc_set.z_max, :z_max)
    end
    
    # 2. Inlets / Outlets (Overwrite External potentially)
    # TODO: Implement inlet logic
    
    # 3. Internal Boundaries
    for ib in bc_set.internal_boundaries
         # Implemented logic for Internal BCs would go here
    end
end

function apply_boundary_conditions!(
    buffers::CFDBuffers,
    grid::GridData,
    bc_set::BoundaryConditionSet,
    dt::Float64,
    par::String
)
    # Apply Velocity BCs
    apply_velocity_bcs!(buffers.u, buffers.v, buffers.w, grid, bc_set, dt)
    
    # Apply Periodic Pressure BCs
    apply_periodic_pressure!(buffers.p, grid, bc_set)
end

"""
    apply_periodic_pressure!(p, grid, bc_set)

BoundaryConditionSetに基づいて圧力の周期境界条件を適用する。
"""
function apply_periodic_pressure!(
    p::Array{Float64, 3},
    grid::GridData,
    bc_set::BoundaryConditionSet
)
    # Y-direction periodic
    if bc_set.y_min.velocity_type == Periodic && bc_set.y_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :y)
    end
    
    # X-direction periodic
    if bc_set.x_min.velocity_type == Periodic && bc_set.x_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :x)
    end
    
    # Z-direction periodic
    if bc_set.z_min.velocity_type == Periodic && bc_set.z_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :z)
    end
end

end # module BoundaryConditions
