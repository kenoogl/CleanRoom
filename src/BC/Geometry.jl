module Geometry

using JSON3
using ..Common
using ..Grid

export GeometryObject, load_geometry, fill_mask!, apply_object_velocity!

"""
    GeometryObject

物体形状定義。
"""
struct GeometryObject
    name::String
    type::Symbol          # :box, :cylinder, :sphere
    params::Dict{Symbol, Any}
    velocity::NTuple{3, Float64}
end

"""
    load_geometry(filepath::String, dim_params::DimensionParams)::Vector{GeometryObject}

Geometry JSONから物体リストを読み込む。
"""
function load_geometry(filepath::String, dim_params::DimensionParams)::Vector{GeometryObject}
    if !isfile(filepath)
         if filepath == "" 
             return GeometryObject[]
         end
         error("Geometry file not found: $filepath")
    end

    L0 = dim_params.L0
    U0 = dim_params.U0

    json_str = read(filepath, String)
    data = JSON3.read(json_str)
    if !haskey(data, :objects)
        error("Geometry file must contain \"objects\".")
    end

    objects = GeometryObject[]
    for item in data.objects
        name = String(item.name)
        type = Symbol(item.type)
        vel_dim = haskey(item, :velocity) ? (Float64(item.velocity[1]), Float64(item.velocity[2]), Float64(item.velocity[3])) : (0.0, 0.0, 0.0)
        vel = vel_dim ./ U0
        
        # Parse params based on type
        params = Dict{Symbol, Any}()
        if type == :box
            pmin = (Float64(item.min[1]), Float64(item.min[2]), Float64(item.min[3]))
            pmax = (Float64(item.max[1]), Float64(item.max[2]), Float64(item.max[3]))
            params[:min] = pmin ./ L0
            params[:max] = pmax ./ L0
        elseif type == :cylinder
            cent = (Float64(item.center[1]), Float64(item.center[2]), Float64(item.center[3]))
            params[:center] = cent ./ L0
            params[:radius] = Float64(item.radius) / L0
            params[:height] = Float64(item.height) / L0
            params[:axis] = Symbol(item.axis)
        elseif type == :sphere
            cent = (Float64(item.center[1]), Float64(item.center[2]), Float64(item.center[3]))
            params[:center] = cent ./ L0
            params[:radius] = Float64(item.radius) / L0
        else
            error("Unknown geometry type: $(type)")
        end

        push!(objects, GeometryObject(name, type, params, vel))
    end
    
    return objects
end

"""
    in_box(x, y, z, min_p, max_p)

直方体内外判定。
"""
function in_box(x, y, z, min_p, max_p)
    return (min_p[1] <= x <= max_p[1]) &&
           (min_p[2] <= y <= max_p[2]) &&
           (min_p[3] <= z <= max_p[3])
end

"""
    in_cylinder(x, y, z, center, radius, height, axis)

円筒内外判定。
"""
function in_cylinder(x, y, z, center, radius, height, axis)
    dx, dy, dz = x - center[1], y - center[2], z - center[3]
    r2 = 0.0
    h_ok = false
    
    if axis == :z
        r2 = dx^2 + dy^2
        h_ok = 0 <= dz <= height
        # Support center as bottom center? "center" usually implies center of bottom or center of volume?
        # Standard: bottom center.
        # But if centered, center - h/2 to center + h/2?
        # Let's assume bottom center for now as it's common in simplifications.
        # Or check requirements... "cylindrical" usually specifies this.
        # Design doc doesn't specify. I'll assume center is min-z (start of cylinder).
    elseif axis == :y
        r2 = dx^2 + dz^2
        h_ok = 0 <= dy <= height
    elseif axis == :x
        r2 = dy^2 + dz^2
        h_ok = 0 <= dx <= height
    end
    
    return h_ok && (r2 <= radius^2)
end

"""
    in_sphere(x, y, z, center, radius)

球内外判定。
"""
function in_sphere(x, y, z, center, radius)
    dist2 = (x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2
    return dist2 <= radius^2
end

function is_inside_object(obj::GeometryObject, x, y, z)
    if obj.type == :box
        return in_box(x, y, z, obj.params[:min], obj.params[:max])
    elseif obj.type == :cylinder
        return in_cylinder(x, y, z, obj.params[:center], obj.params[:radius], obj.params[:height], obj.params[:axis])
    elseif obj.type == :sphere
        return in_sphere(x, y, z, obj.params[:center], obj.params[:radius])
    end
    return false
end


"""
    fill_mask!(mask, objects, grid, par)

マスク配列を生成する（流体=1, 物体=0）。
"""
function fill_mask!(
    mask::Array{Float64, 3},
    objects::Vector{GeometryObject},
    grid::GridData,
    par::String # parallel param, unused in this sequential loop for now
)
    # Reset mask to 1.0 (Fluid)
    fill!(mask, 1.0)
    
    mx, my, mz = grid.mx, grid.my, grid.mz
    
    # Iterate over all cells (including ghosts? Usually ghosts handle BCs, but mask applies to domain)
    # Mask should cover ghost cells if object extends there.
    
    for k in 1:mz
        z = grid.z_center[k]
        for j in 1:my
            y = grid.y[j]
            for i in 1:mx
                x = grid.x[i]
                
                is_solid = false
                for obj in objects
                    if is_inside_object(obj, x, y, z)
                        is_solid = true
                        break
                    end
                end
                
                if is_solid
                    mask[i, j, k] = 0.0
                end
            end
        end
    end
end

function apply_object_velocity!(
    u::Array{Float64, 3},
    v::Array{Float64, 3},
    w::Array{Float64, 3},
    objects::Vector{GeometryObject},
    grid::GridData
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    @inbounds for k in 1:mz
        z = grid.z_center[k]
        for j in 1:my
            y = grid.y[j]
            for i in 1:mx
                x = grid.x[i]
                for obj in objects
                    if is_inside_object(obj, x, y, z)
                        u[i, j, k] = obj.velocity[1]
                        v[i, j, k] = obj.velocity[2]
                        w[i, j, k] = obj.velocity[3]
                        break
                    end
                end
            end
        end
    end
end

end # module Geometry
