module BoundaryConditions

using ..Common
using ..Grid
using ..Fields

export VelocityBCType, ExternalBC, InletOutlet, InternalBoundary, BoundaryConditionSet
export Dirichlet, Neumann, Outflow, Periodic, Symmetric, Wall, SlidingWall, Inlet 
export apply_boundary_conditions!, apply_outflow!, apply_velocity_bcs!
export apply_periodic_velocity!, apply_periodic_pressure!, apply_periodic_face_velocity!
export apply_face_velocity_bcs!, apply_boundary_mask!, apply_pressure_bcs!
export apply_outflow_region!, update_outflow_mask!

@enum VelocityBCType begin
    Inlet
    Dirichlet
    Neumann
    Outflow
    Periodic
    Symmetric
    Wall
    SlidingWall
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

function compute_region_average_velocity(
    u::Array{Float64, 3},
    v::Array{Float64, 3},
    w::Array{Float64, 3},
    grid::GridData,
    face::Symbol,
    region_check::Function
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    sum_vel = 0.0
    count = 0
    if face == :x_min || face == :x_max
        i = (face == :x_min) ? 3 : mx-2
        for k in 3:mz-2, j in 3:my-2
            if region_check(i, j, k)
                sum_vel += u[i, j, k]
                count += 1
            end
        end
    elseif face == :y_min || face == :y_max
        j = (face == :y_min) ? 3 : my-2
        for k in 3:mz-2, i in 3:mx-2
            if region_check(i, j, k)
                sum_vel += v[i, j, k]
                count += 1
            end
        end
    elseif face == :z_min || face == :z_max
        k = (face == :z_min) ? 3 : mz-2
        for j in 3:my-2, i in 3:mx-2
            if region_check(i, j, k)
                sum_vel += w[i, j, k]
                count += 1
            end
        end
    end
    return (count > 0) ? sum_vel / count : 0.0
end

@inline function face_from_normal(normal::NTuple{3, Int})::Symbol
    if normal == (1, 0, 0)
        return :x_max
    elseif normal == (-1, 0, 0)
        return :x_min
    elseif normal == (0, 1, 0)
        return :y_max
    elseif normal == (0, -1, 0)
        return :y_min
    elseif normal == (0, 0, 1)
        return :z_max
    elseif normal == (0, 0, -1)
        return :z_min
    end
    error("normal must be axis-aligned: $(normal)")
end

@inline function is_in_rectangular_patch(
    x::Float64,
    y::Float64,
    z::Float64,
    position::NTuple{3, Float64},
    size::NTuple{2, Float64},
    normal::NTuple{3, Int},
    grid::GridData,
    k::Int
)::Bool
    if normal[1] != 0
        return abs(x - position[1]) <= grid.dx &&
               abs(y - position[2]) <= size[1] * 0.5 &&
               abs(z - position[3]) <= size[2] * 0.5
    elseif normal[2] != 0
        return abs(y - position[2]) <= grid.dy &&
               abs(x - position[1]) <= size[1] * 0.5 &&
               abs(z - position[3]) <= size[2] * 0.5
    elseif normal[3] != 0
        return abs(z - position[3]) <= grid.dz[k] &&
               abs(x - position[1]) <= size[1] * 0.5 &&
               abs(y - position[2]) <= size[2] * 0.5
    end
    return false
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
    apply_outflow_region!(phi, grid, face, dt, Uc, region_check)

領域指定付きの対流流出条件を適用する。
`Outflow` 条件が指定された `InletOutlet`（パッチ）に対して、
`region_check` で指定された範囲内のゴーストセルのみを更新する。

- `Uc`: 領域内の平均流出速度（`compute_region_average_velocity` で計算）
- `region_check(i, j, k)`: インデックスがパッチ内かどうかを判定する関数
"""
function apply_outflow_region!(
    phi::Array{Float64, 3},
    grid::GridData,
    face::Symbol,
    dt::Float64,
    Uc::Float64,
    region_check::Function
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    if face == :x_max
        dx = grid.dx
        c = Uc * dt / dx
        i_inner = mx - 2
        @inbounds for k in 3:mz-2, j in 3:my-2
            if region_check(i_inner, j, k)
                phi[mx-1, j, k] = phi[mx-1, j, k] - c * (phi[mx-1, j, k] - phi[mx-2, j, k])
                phi[mx, j, k] = phi[mx-1, j, k]
            end
        end
    elseif face == :x_min
        dx = grid.dx
        c = Uc * dt / dx
        i_inner = 3
        @inbounds for k in 3:mz-2, j in 3:my-2
            if region_check(i_inner, j, k)
                phi[2, j, k] = phi[2, j, k] - c * (phi[3, j, k] - phi[2, j, k])
                phi[1, j, k] = phi[2, j, k]
            end
        end
    elseif face == :y_max
        dy = grid.dy
        c = Uc * dt / dy
        j_inner = my - 2
        @inbounds for k in 3:mz-2, i in 3:mx-2
            if region_check(i, j_inner, k)
                phi[i, my-1, k] = phi[i, my-1, k] - c * (phi[i, my-1, k] - phi[i, my-2, k])
                phi[i, my, k] = phi[i, my-1, k]
            end
        end
    elseif face == :y_min
        dy = grid.dy
        c = Uc * dt / dy
        j_inner = 3
        @inbounds for k in 3:mz-2, i in 3:mx-2
            if region_check(i, j_inner, k)
                phi[i, 2, k] = phi[i, 2, k] - c * (phi[i, 3, k] - phi[i, 2, k])
                phi[i, 1, k] = phi[i, 2, k]
            end
        end
    elseif face == :z_max
        dz = grid.dz[mz-2]
        c = Uc * dt / dz
        k_inner = mz - 2
        @inbounds for j in 3:my-2, i in 3:mx-2
            if region_check(i, j, k_inner)
                phi[i, j, mz-1] = phi[i, j, mz-1] - c * (phi[i, j, mz-1] - phi[i, j, mz-2])
                phi[i, j, mz] = phi[i, j, mz-1]
            end
        end
    elseif face == :z_min
        dz = grid.dz[3]
        c = Uc * dt / dz
        k_inner = 3
        @inbounds for j in 3:my-2, i in 3:mx-2
            if region_check(i, j, k_inner)
                phi[i, j, 2] = phi[i, j, 2] - c * (phi[i, j, 3] - phi[i, j, 2])
                phi[i, j, 1] = phi[i, j, 2]
            end
        end
    end
end

"""
    apply_symmetric_bc!(u, v, w, grid, face)

対称境界条件を適用する。
- 法線方向速度成分: 符号反転（u_n = 0）
- 接線方向速度成分: ミラーリング（∂u_t/∂n = 0）
"""
function apply_symmetric_bc!(
    u::Array{Float64, 3},
    v::Array{Float64, 3},
    w::Array{Float64, 3},
    grid::GridData,
    face::Symbol
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    if face == :x_min
        # 境界: i=3, ゴースト: i=1,2
        # 法線方向: u (符号反転)
        # 接線方向: v, w (ミラーリング)
        @inbounds for k in 1:mz, j in 1:my
            u[2, j, k] = -u[3, j, k]
            u[1, j, k] = -u[4, j, k]
            v[2, j, k] = v[3, j, k]
            v[1, j, k] = v[4, j, k]
            w[2, j, k] = w[3, j, k]
            w[1, j, k] = w[4, j, k]
        end
    elseif face == :x_max
        # 境界: i=mx-2, ゴースト: i=mx-1,mx
        @inbounds for k in 1:mz, j in 1:my
            u[mx-1, j, k] = -u[mx-2, j, k]
            u[mx, j, k] = -u[mx-3, j, k]
            v[mx-1, j, k] = v[mx-2, j, k]
            v[mx, j, k] = v[mx-3, j, k]
            w[mx-1, j, k] = w[mx-2, j, k]
            w[mx, j, k] = w[mx-3, j, k]
        end
    elseif face == :y_min
        # 境界: j=3, ゴースト: j=1,2
        # 法線方向: v (符号反転)
        # 接線方向: u, w (ミラーリング)
        @inbounds for k in 1:mz, i in 1:mx
            v[i, 2, k] = -v[i, 3, k]
            v[i, 1, k] = -v[i, 4, k]
            u[i, 2, k] = u[i, 3, k]
            u[i, 1, k] = u[i, 4, k]
            w[i, 2, k] = w[i, 3, k]
            w[i, 1, k] = w[i, 4, k]
        end
    elseif face == :y_max
        # 境界: j=my-2, ゴースト: j=my-1,my
        @inbounds for k in 1:mz, i in 1:mx
            v[i, my-1, k] = -v[i, my-2, k]
            v[i, my, k] = -v[i, my-3, k]
            u[i, my-1, k] = u[i, my-2, k]
            u[i, my, k] = u[i, my-3, k]
            w[i, my-1, k] = w[i, my-2, k]
            w[i, my, k] = w[i, my-3, k]
        end
    elseif face == :z_min
        # 境界: k=3, ゴースト: k=1,2
        # 法線方向: w (符号反転)
        # 接線方向: u, v (ミラーリング)
        @inbounds for j in 1:my, i in 1:mx
            w[i, j, 2] = -w[i, j, 3]
            w[i, j, 1] = -w[i, j, 4]
            u[i, j, 2] = u[i, j, 3]
            u[i, j, 1] = u[i, j, 4]
            v[i, j, 2] = v[i, j, 3]
            v[i, j, 1] = v[i, j, 4]
        end
    elseif face == :z_max
        # 境界: k=mz-2, ゴースト: k=mz-1,mz
        @inbounds for j in 1:my, i in 1:mx
            w[i, j, mz-1] = -w[i, j, mz-2]
            w[i, j, mz] = -w[i, j, mz-3]
            u[i, j, mz-1] = u[i, j, mz-2]
            u[i, j, mz] = u[i, j, mz-3]
            v[i, j, mz-1] = v[i, j, mz-2]
            v[i, j, mz] = v[i, j, mz-3]
        end
    end
end

@inline function in_internal_cylinder(
    x::Float64,
    y::Float64,
    z::Float64,
    center::NTuple{3, Float64},
    radius::Float64,
    height::Float64,
    axis::Symbol
)::Bool
    dx, dy, dz = x - center[1], y - center[2], z - center[3]
    r2 = 0.0
    h_ok = false
    if axis == :z
        r2 = dx^2 + dy^2
        h_ok = 0.0 <= dz <= height
    elseif axis == :y
        r2 = dx^2 + dz^2
        h_ok = 0.0 <= dy <= height
    elseif axis == :x
        r2 = dy^2 + dz^2
        h_ok = 0.0 <= dx <= height
    end
    return h_ok && (r2 <= radius^2)
end

@inline function is_in_internal_boundary(
    ib::InternalBoundary,
    x::Float64,
    y::Float64,
    z::Float64
)::Bool
    if ib.type == :rectangular
        return (ib.region_min[1] <= x <= ib.region_max[1]) &&
               (ib.region_min[2] <= y <= ib.region_max[2]) &&
               (ib.region_min[3] <= z <= ib.region_max[3])
    elseif ib.type == :cylindrical
        return in_internal_cylinder(x, y, z, ib.center, ib.radius, ib.height, ib.axis)
    end
    return false
end

"""
    set_boundary_value!(phi, grid, mask, face, val, type)

Helper to set Dirichlet/Neumann on a face, considering the mask.
For Dirichlet flow, the velocity is masked to 0 on solid regions.
"""
function set_boundary_value!(
    phi::Array{Float64, 3},
    grid::GridData,
    mask::Array{Float64, 3},
    face::Symbol,
    val::Float64,
    type::VelocityBCType
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    if type == Dirichlet || type == Inlet
        if face == :x_min
             # Use mask from first internal cell (i=3)
             @inbounds for k in 1:mz, j in 1:my
                 phi[1, j, k] = val * mask[3, j, k]
                 phi[2, j, k] = val * mask[3, j, k]
             end
        elseif face == :x_max
             @inbounds for k in 1:mz, j in 1:my
                 phi[mx-1, j, k] = val * mask[mx-2, j, k]
                 phi[mx, j, k] = val * mask[mx-2, j, k]
             end
        elseif face == :y_min
             @inbounds for k in 1:mz, i in 1:mx
                 phi[i, 1, k] = val * mask[i, 3, k]
                 phi[i, 2, k] = val * mask[i, 3, k]
             end
        elseif face == :y_max
             @inbounds for k in 1:mz, i in 1:mx
                 phi[i, my-1, k] = val * mask[i, my-2, k]
                 phi[i, my, k] = val * mask[i, my-2, k]
             end
        elseif face == :z_min
             @inbounds for j in 1:my, i in 1:mx
                 phi[i, j, 1] = val * mask[i, j, 3]
                 phi[i, j, 2] = val * mask[i, j, 3]
             end
        elseif face == :z_max
             @inbounds for j in 1:my, i in 1:mx
                 phi[i, j, mz-1] = val * mask[i, j, mz-2]
                 phi[i, j, mz] = val * mask[i, j, mz-2]
             end
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
    apply_velocity_bcs!(u, v, w, grid, mask, bc_set, dt)

指定された速度場(u, v, w)に対して境界条件を適用する。
"""
function apply_velocity_bcs!(
    u::Array{Float64, 3},
    v::Array{Float64, 3},
    w::Array{Float64, 3},
    grid::GridData,
    mask::Array{Float64, 3},
    bc_set::BoundaryConditionSet,
    dt::Float64
)
    # 1. External Boundaries
    function apply_ext(bc, face)
        if bc.velocity_type == Dirichlet || bc.velocity_type == Inlet
             set_boundary_value!(u, grid, mask, face, bc.velocity_value[1], bc.velocity_type)
             set_boundary_value!(v, grid, mask, face, bc.velocity_value[2], bc.velocity_type)
             set_boundary_value!(w, grid, mask, face, bc.velocity_value[3], bc.velocity_type)
        elseif bc.velocity_type == Neumann
             set_boundary_value!(u, grid, mask, face, 0.0, Neumann)
             set_boundary_value!(v, grid, mask, face, 0.0, Neumann)
             set_boundary_value!(w, grid, mask, face, 0.0, Neumann)
        elseif bc.velocity_type == Outflow
             Uc = compute_face_average_velocity(u, v, w, grid, face)
             apply_outflow!(u, grid, face, dt, Uc)
             apply_outflow!(v, grid, face, dt, Uc)
             apply_outflow!(w, grid, face, dt, Uc)
        elseif bc.velocity_type == Symmetric
             apply_symmetric_bc!(u, v, w, grid, face)
        elseif bc.velocity_type == Wall
             set_boundary_value!(u, grid, mask, face, 0.0, Dirichlet)
             set_boundary_value!(v, grid, mask, face, 0.0, Dirichlet)
             set_boundary_value!(w, grid, mask, face, 0.0, Dirichlet)
        elseif bc.velocity_type == SlidingWall
             set_boundary_value!(u, grid, mask, face, bc.velocity_value[1], Dirichlet)
             set_boundary_value!(v, grid, mask, face, bc.velocity_value[2], Dirichlet)
             set_boundary_value!(w, grid, mask, face, bc.velocity_value[3], Dirichlet)
        end
    end

    # Check for periodic BC pairs
    y_periodic = bc_set.y_min.velocity_type == Periodic && bc_set.y_max.velocity_type == Periodic
    x_periodic = bc_set.x_min.velocity_type == Periodic && bc_set.x_max.velocity_type == Periodic
    z_periodic = bc_set.z_min.velocity_type == Periodic && bc_set.z_max.velocity_type == Periodic

    if (bc_set.y_min.velocity_type == Periodic) != (bc_set.y_max.velocity_type == Periodic)
        error("Periodic boundary requires both y_min and y_max to be periodic.")
    end
    if (bc_set.x_min.velocity_type == Periodic) != (bc_set.x_max.velocity_type == Periodic)
        error("Periodic boundary requires both x_min and x_max to be periodic.")
    end
    if (bc_set.z_min.velocity_type == Periodic) != (bc_set.z_max.velocity_type == Periodic)
        error("Periodic boundary requires both z_min and z_max to be periodic.")
    end
    
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
    
    # 2. Inlets
    for inlet in bc_set.inlets
        for k in 3:grid.mz-2, j in 3:grid.my-2, i in 3:grid.mx-2
            x, y, z = grid.x[i], grid.y[j], grid.z_center[k]
            if inlet.type != :rectangular
                error("Unsupported inlet type: $(inlet.type)")
            end
            in_range = is_in_rectangular_patch(
                x, y, z,
                inlet.position,
                inlet.size,
                inlet.normal,
                grid,
                k
            )
            if in_range
                u[i, j, k] = inlet.velocity[1]
                v[i, j, k] = inlet.velocity[2]
                w[i, j, k] = inlet.velocity[3]
            end
        end
    end

    # 3. Outlets
    for outlet in bc_set.outlets
        if outlet.type != :rectangular
            error("Unsupported outlet type: $(outlet.type)")
        end
        face = face_from_normal(outlet.normal)
        region_check = (i, j, k) -> is_in_rectangular_patch(
            grid.x[i],
            grid.y[j],
            grid.z_center[k],
            outlet.position,
            outlet.size,
            outlet.normal,
            grid,
            k
        )
        if outlet.condition == Dirichlet
            for k in 3:grid.mz-2, j in 3:grid.my-2, i in 3:grid.mx-2
                if region_check(i, j, k)
                    u[i, j, k] = outlet.velocity[1]
                    v[i, j, k] = outlet.velocity[2]
                    w[i, j, k] = outlet.velocity[3]
                end
            end
        elseif outlet.condition == Outflow
            Uc = compute_region_average_velocity(u, v, w, grid, face, region_check)
            apply_outflow_region!(u, grid, face, dt, Uc, region_check)
            apply_outflow_region!(v, grid, face, dt, Uc, region_check)
            apply_outflow_region!(w, grid, face, dt, Uc, region_check)
        end
    end
    
    # 4. Internal Boundaries
    for ib in bc_set.internal_boundaries
        for k in 3:grid.mz-2, j in 3:grid.my-2, i in 3:grid.mx-2
            x, y, z = grid.x[i], grid.y[j], grid.z_center[k]
            if is_in_internal_boundary(ib, x, y, z)
                u[i, j, k] = ib.velocity[1]
                v[i, j, k] = ib.velocity[2]
                w[i, j, k] = ib.velocity[3]
            end
        end
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
    apply_velocity_bcs!(buffers.u, buffers.v, buffers.w, grid, buffers.mask, bc_set, dt)

    # Apply Pressure BCs (Neumann / Periodic)
    apply_pressure_bcs!(buffers.p, grid, buffers.mask, bc_set)
    apply_internal_pressure_bcs!(buffers.p, grid, bc_set)
end

function apply_pressure_bcs!(
    p::Array{Float64, 3},
    grid::GridData,
    mask::Array{Float64, 3},
    bc_set::BoundaryConditionSet
)
    # Y-direction
    if bc_set.y_min.velocity_type == Periodic && bc_set.y_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :y)
    else
        # y_min
        p_type_min = (bc_set.y_min.velocity_type == Outflow) ? Dirichlet : Neumann
        set_boundary_value!(p, grid, mask, :y_min, 0.0, p_type_min)
        # y_max
        p_type_max = (bc_set.y_max.velocity_type == Outflow) ? Dirichlet : Neumann
        set_boundary_value!(p, grid, mask, :y_max, 0.0, p_type_max)
    end
    # X-direction
    if bc_set.x_min.velocity_type == Periodic && bc_set.x_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :x)
    else
        # x_min
        p_type_min = (bc_set.x_min.velocity_type == Outflow) ? Dirichlet : Neumann
        set_boundary_value!(p, grid, mask, :x_min, 0.0, p_type_min)
        # x_max
        p_type_max = (bc_set.x_max.velocity_type == Outflow) ? Dirichlet : Neumann
        set_boundary_value!(p, grid, mask, :x_max, 0.0, p_type_max)
    end
    # Z-direction
    if bc_set.z_min.velocity_type == Periodic && bc_set.z_max.velocity_type == Periodic
        apply_periodic_pressure!(p, grid, :z)
    else
        # z_min
        p_type_min = (bc_set.z_min.velocity_type == Outflow) ? Dirichlet : Neumann
        set_boundary_value!(p, grid, mask, :z_min, 0.0, p_type_min)
        # z_max
        p_type_max = (bc_set.z_max.velocity_type == Outflow) ? Dirichlet : Neumann
        set_boundary_value!(p, grid, mask, :z_max, 0.0, p_type_max)
    end
end

function apply_internal_pressure_bcs!(
    p::Array{Float64, 3},
    grid::GridData,
    bc_set::BoundaryConditionSet
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    for ib in bc_set.internal_boundaries
        ni, nj, nk = ib.normal
        di = (ni != 0) ? -ni : 0
        dj = (nj != 0) ? -nj : 0
        dk = (nk != 0) ? -nk : 0
        @inbounds for k in 3:mz-2, j in 3:my-2, i in 3:mx-2
            x, y, z = grid.x[i], grid.y[j], grid.z_center[k]
            if is_in_internal_boundary(ib, x, y, z)
                src_i = i + di
                src_j = j + dj
                src_k = k + dk
                if 1 <= src_i <= mx && 1 <= src_j <= my && 1 <= src_k <= mz
                    p[i, j, k] = p[src_i, src_j, src_k]
                end
            end
        end
    end
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

"""
    apply_periodic_face_velocity!(u_face, v_face, w_face, grid, axis)

セルフェイス速度の周期境界条件を適用する。
セルフェイス速度はインデックスがセルセンターと異なるため専用の処理が必要。
"""
function apply_periodic_face_velocity!(
    u_face::Array{Float64, 3},
    v_face::Array{Float64, 3},
    w_face::Array{Float64, 3},
    grid::GridData,
    axis::Symbol
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    if axis == :y
        # Y方向周期境界
        # u_face_x, w_face_z: セルセンターと同じインデックス
        @inbounds for k in 1:mz, i in 1:mx
            u_face[i, 1, k] = u_face[i, my-3, k]
            u_face[i, 2, k] = u_face[i, my-2, k]
            u_face[i, my-1, k] = u_face[i, 3, k]
            u_face[i, my, k] = u_face[i, 4, k]
            
            w_face[i, 1, k] = w_face[i, my-3, k]
            w_face[i, 2, k] = w_face[i, my-2, k]
            w_face[i, my-1, k] = w_face[i, 3, k]
            w_face[i, my, k] = w_face[i, 4, k]
        end
        
        # v_face_y: j+1/2位置なのでインデックスがシフト
        @inbounds for k in 1:mz, i in 1:mx
            v_face[i, 1, k] = v_face[i, my-3, k]
            v_face[i, 2, k] = v_face[i, my-2, k]
            v_face[i, my-1, k] = v_face[i, 3, k]
            v_face[i, my, k] = v_face[i, 4, k]
        end
        
    elseif axis == :x
        # X方向周期境界
        # v_face_y, w_face_z: セルセンターと同じインデックス
        @inbounds for k in 1:mz, j in 1:my
            v_face[1, j, k] = v_face[mx-3, j, k]
            v_face[2, j, k] = v_face[mx-2, j, k]
            v_face[mx-1, j, k] = v_face[3, j, k]
            v_face[mx, j, k] = v_face[4, j, k]
            
            w_face[1, j, k] = w_face[mx-3, j, k]
            w_face[2, j, k] = w_face[mx-2, j, k]
            w_face[mx-1, j, k] = w_face[3, j, k]
            w_face[mx, j, k] = w_face[4, j, k]
        end
        
        # u_face_x: i+1/2位置なのでインデックスがシフト
        @inbounds for k in 1:mz, j in 1:my
            u_face[1, j, k] = u_face[mx-3, j, k]
            u_face[2, j, k] = u_face[mx-2, j, k]
            u_face[mx-1, j, k] = u_face[3, j, k]
            u_face[mx, j, k] = u_face[4, j, k]
        end
        
    elseif axis == :z
        # Z方向周期境界
        # u_face_x, v_face_y: セルセンターと同じインデックス
        @inbounds for j in 1:my, i in 1:mx
            u_face[i, j, 1] = u_face[i, j, mz-3]
            u_face[i, j, 2] = u_face[i, j, mz-2]
            u_face[i, j, mz-1] = u_face[i, j, 3]
            u_face[i, j, mz] = u_face[i, j, 4]
            
            v_face[i, j, 1] = v_face[i, j, mz-3]
            v_face[i, j, 2] = v_face[i, j, mz-2]
            v_face[i, j, mz-1] = v_face[i, j, 3]
            v_face[i, j, mz] = v_face[i, j, 4]
        end
        
        # w_face_z: k+1/2位置なのでインデックスがシフト
        @inbounds for j in 1:my, i in 1:mx
            w_face[i, j, 1] = w_face[i, j, mz-3]
            w_face[i, j, 2] = w_face[i, j, mz-2]
            w_face[i, j, mz-1] = w_face[i, j, 3]
            w_face[i, j, mz] = w_face[i, j, 4]
        end
    end
end


"""
    apply_symmetric_face_bc!(u_f, v_f, w_f, grid, face)

セルフェイス速度に対称境界条件を適用する。
- 法線方向フェイス速度: 境界位置で 0、ゴーストで符号反転
- 接線方向フェイス速度: ゴーストでミラーリング
"""
function apply_symmetric_face_bc!(
    u_f::Array{Float64, 3},
    v_f::Array{Float64, 3},
    w_f::Array{Float64, 3},
    grid::GridData,
    face::Symbol
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    if face == :x_min
        # 境界位置: i=3 (x=2.5)
        @inbounds for k in 1:mz, j in 1:my
            u_f[3, j, k] = 0.0
            u_f[2, j, k] = -u_f[4, j, k]
            u_f[1, j, k] = -u_f[5, j, k]
            v_f[2, j, k] = v_f[3, j, k]
            v_f[1, j, k] = v_f[4, j, k]
            w_f[2, j, k] = w_f[3, j, k]
            w_f[1, j, k] = w_f[4, j, k]
        end
    elseif face == :x_max
        # 境界位置: i=mx-1 (x=mx-1.5)
        @inbounds for k in 1:mz, j in 1:my
            u_f[mx-1, j, k] = 0.0
            u_f[mx, j, k] = -u_f[mx-2, j, k]
            v_f[mx-1, j, k] = v_f[mx-2, j, k]
            v_f[mx, j, k] = v_f[mx-3, j, k]
            w_f[mx-1, j, k] = w_f[mx-2, j, k]
            w_f[mx, j, k] = w_f[mx-3, j, k]
        end
    elseif face == :y_min
        # 境界位置: j=3 (y=2.5)
        @inbounds for k in 1:mz, i in 1:mx
            v_f[i, 3, k] = 0.0
            v_f[i, 2, k] = -v_f[i, 4, k]
            v_f[i, 1, k] = -v_f[i, 5, k]
            u_f[i, 2, k] = u_f[i, 3, k]
            u_f[i, 1, k] = u_f[i, 4, k]
            w_f[i, 2, k] = w_f[i, 3, k]
            w_f[i, 1, k] = w_f[i, 4, k]
        end
    elseif face == :y_max
        # 境界位置: j=my-1 (y=my-1.5)
        @inbounds for k in 1:mz, i in 1:mx
            v_f[i, my-1, k] = 0.0
            v_f[i, my, k] = -v_f[i, my-2, k]
            u_f[i, my-1, k] = u_f[i, my-2, k]
            u_f[i, my, k] = u_f[i, my-3, k]
            w_f[i, my-1, k] = w_f[i, my-2, k]
            w_f[i, my, k] = w_f[i, my-3, k]
        end
    elseif face == :z_min
        # 境界位置: k=3 (z=2.5)
        @inbounds for j in 1:my, i in 1:mx
            w_f[i, j, 3] = 0.0
            w_f[i, j, 2] = -w_f[i, j, 4]
            w_f[i, j, 1] = -w_f[i, j, 5]
            u_f[i, j, 2] = u_f[i, j, 3]
            u_f[i, j, 1] = u_f[i, j, 4]
            v_f[i, j, 2] = v_f[i, j, 3]
            v_f[i, j, 1] = v_f[i, j, 4]
        end
    elseif face == :z_max
        # 境界位置: k=mz-1 (z=mz-1.5)
        @inbounds for j in 1:my, i in 1:mx
            w_f[i, j, mz-1] = 0.0
            w_f[i, j, mz] = -w_f[i, j, mz-2]
            u_f[i, j, mz-1] = u_f[i, j, mz-2]
            u_f[i, j, mz] = u_f[i, j, mz-3]
            v_f[i, j, mz-1] = v_f[i, j, mz-2]
            v_f[i, j, mz] = v_f[i, j, mz-3]
        end
    end
end

"""
    apply_face_velocity_bcs!(u_f, v_f, w_f, grid, mask, bc_set, dt)

セルフェイス速度の境界条件を一括適用する。
"""
function apply_face_velocity_bcs!(
    u_f::Array{Float64, 3},
    v_f::Array{Float64, 3},
    w_f::Array{Float64, 3},
    grid::GridData,
    mask::Array{Float64, 3},
    bc_set::BoundaryConditionSet,
    dt::Float64
)
    # 1. External Boundaries
    function apply_ext_face(bc, face)
        if bc.velocity_type == Symmetric
            apply_symmetric_face_bc!(u_f, v_f, w_f, grid, face)
        elseif bc.velocity_type == Wall
            # Handled by mask in FractionalStep, but we can be explicit here
            # for boundary faces.
            if face == :x_min; u_f[3, :, :] .= 0.0; elseif face == :x_max; u_f[grid.mx-1, :, :] .= 0.0; end
            if face == :y_min; v_f[:, 3, :] .= 0.0; elseif face == :y_max; v_f[:, grid.my-1, :] .= 0.0; end
            if face == :z_min; w_f[:, :, 3] .= 0.0; elseif face == :z_max; w_f[:, :, grid.mz-1] .= 0.0; end
        elseif bc.velocity_type == SlidingWall
            if face == :x_min; u_f[3, :, :] .= bc.velocity_value[1]; elseif face == :x_max; u_f[grid.mx-1, :, :] .= bc.velocity_value[1]; end
            if face == :y_min; v_f[:, 3, :] .= bc.velocity_value[2]; elseif face == :y_max; v_f[:, grid.my-1, :] .= bc.velocity_value[2]; end
            if face == :z_min; w_f[:, :, 3] .= bc.velocity_value[3]; elseif face == :z_max; w_f[:, :, grid.mz-1] .= bc.velocity_value[3]; end
        end
    end

    x_periodic = bc_set.x_min.velocity_type == Periodic && bc_set.x_max.velocity_type == Periodic
    y_periodic = bc_set.y_min.velocity_type == Periodic && bc_set.y_max.velocity_type == Periodic
    z_periodic = bc_set.z_min.velocity_type == Periodic && bc_set.z_max.velocity_type == Periodic

    if x_periodic
        apply_periodic_face_velocity!(u_f, v_f, w_f, grid, :x)
    else
        apply_ext_face(bc_set.x_min, :x_min)
        apply_ext_face(bc_set.x_max, :x_max)
    end

    if y_periodic
        apply_periodic_face_velocity!(u_f, v_f, w_f, grid, :y)
    else
        apply_ext_face(bc_set.y_min, :y_min)
        apply_ext_face(bc_set.y_max, :y_max)
    end

    if z_periodic
        apply_periodic_face_velocity!(u_f, v_f, w_f, grid, :z)
    else
        apply_ext_face(bc_set.z_min, :z_min)
        apply_ext_face(bc_set.z_max, :z_max)
    end
end


"""
    apply_boundary_mask!(mask, grid, bc_set)

外部境界条件に基づいてゴーストセルのマスク値を設定する。
- Wall, Symmetric: 0.0 (Solid) -> Flux blocked
- Inlet, Outflow, Neumann: 1.0 (Fluid) -> Flux allowed (handled by ghost cells)
"""
function apply_boundary_mask!(
    mask::Array{Float64, 3},
    grid::GridData,
    bc_set::BoundaryConditionSet
)
    mx, my, mz = grid.mx, grid.my, grid.mz

    function set_mask(bc, face, val)
        if face == :x_min
            mask[1:2, :, :] .= val
        elseif face == :x_max
            mask[mx-1:mx, :, :] .= val
        elseif face == :y_min
            mask[:, 1:2, :] .= val
        elseif face == :y_max
            mask[:, my-1:my, :] .= val
        elseif face == :z_min
            mask[:, :, 1:2] .= val
        elseif face == :z_max
            mask[:, :, mz-1:mz] .= val
        end
    end

    # Default all ghosts to 1.0 (Fluid) unless Wall/Symmetric
    # Actually, iterate through BCs
    for (face, bc) in [
        (:x_min, bc_set.x_min), (:x_max, bc_set.x_max),
        (:y_min, bc_set.y_min), (:y_max, bc_set.y_max),
        (:z_min, bc_set.z_min), (:z_max, bc_set.z_max)
    ]
        if bc.velocity_type == Wall || bc.velocity_type == Symmetric || bc.velocity_type == SlidingWall
            set_mask(bc, face, 0.0)
        else
            set_mask(bc, face, 1.0)
        end
    end
end

"""
    update_outflow_mask!(mask, grid, bc_set, val)

Outflow境界のゴーストセルマスクを指定した値に変更する。
- 圧力計算時: val = 0.0 (Solid) -> Neumann条件
- 対流計算時: val = 1.0 (Fluid) -> 流れを許容
"""
function update_outflow_mask!(
    mask::Array{Float64, 3},
    grid::GridData,
    bc_set::BoundaryConditionSet,
    val::Float64
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    function set_val(face)
        if face == :x_min
            mask[1:2, :, :] .= val
        elseif face == :x_max
            mask[mx-1:mx, :, :] .= val
        elseif face == :y_min
            mask[:, 1:2, :] .= val
        elseif face == :y_max
            mask[:, my-1:my, :] .= val
        elseif face == :z_min
            mask[:, :, 1:2] .= val
        elseif face == :z_max
            mask[:, :, mz-1:mz] .= val
        end
    end

    if bc_set.x_min.velocity_type == Outflow; set_val(:x_min); end
    if bc_set.x_max.velocity_type == Outflow; set_val(:x_max); end
    if bc_set.y_min.velocity_type == Outflow; set_val(:y_min); end
    if bc_set.y_max.velocity_type == Outflow; set_val(:y_max); end
    if bc_set.z_min.velocity_type == Outflow; set_val(:z_min); end
    if bc_set.z_max.velocity_type == Outflow; set_val(:z_max); end
end

end # module BoundaryConditions
