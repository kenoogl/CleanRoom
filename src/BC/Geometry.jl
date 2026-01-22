module Geometry

using JSON3
using ..Common
using ..Grid

export GeometryObject, load_geometry, fill_mask!

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
    load_geometry(filepath::String)::Vector{GeometryObject}

Geometry JSONから物体リストを読み込む。
"""
function load_geometry(filepath::String)::Vector{GeometryObject}
    if !isfile(filepath)
         # If file doesn't exist, return empty list (no objects)
         # Or error? Spec says Requirement 9.3 "Geometry JSON読込".
         # Usually valid path expected. But if no geometry, user might not provide file?
         # Design doc system flow doesn't explicitly handle optionality.
         # But "Objects" usually exist.
         # I'll return empty if file not found, but log warning?
         # Or error. Let's error to be safe, unless explicit empty.
         # But wait, main driver might pass a default path.
         if filepath == "" 
             return GeometryObject[]
         end
         # If path provided but missing, error.
         error("Geometry file not found: $filepath")
    end

    json_str = read(filepath, String)
    data = JSON3.read(json_str) # Expects list or object with list?
    # Usually list of objects.
    # Assume root is array of objects.
    
    objects = GeometryObject[]
    for item in data
        name = String(item.name)
        type = Symbol(item.type)
        vel = haskey(item, :velocity) ? (Float64(item.velocity[1]), Float64(item.velocity[2]), Float64(item.velocity[3])) : (0.0, 0.0, 0.0)
        
        # Parse params based on type
        params = Dict{Symbol, Any}()
        if type == :box
            params[:min] = (Float64(item.params.min[1]), Float64(item.params.min[2]), Float64(item.params.min[3]))
            params[:max] = (Float64(item.params.max[1]), Float64(item.params.max[2]), Float64(item.params.max[3]))
        elseif type == :cylinder
            params[:center] = (Float64(item.params.center[1]), Float64(item.params.center[2]), Float64(item.params.center[3]))
            params[:radius] = Float64(item.params.radius)
            params[:height] = Float64(item.params.height)
            params[:axis] = Symbol(item.params.axis)
        elseif type == :sphere
            params[:center] = (Float64(item.params.center[1]), Float64(item.params.center[2]), Float64(item.params.center[3]))
            params[:radius] = Float64(item.params.radius)
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
                    if obj.type == :box
                        if in_box(x, y, z, obj.params[:min], obj.params[:max])
                            is_solid = true
                            break
                        end
                    elseif obj.type == :cylinder
                        if in_cylinder(x, y, z, obj.params[:center], obj.params[:radius], obj.params[:height], obj.params[:axis])
                            is_solid = true
                            break
                        end
                    elseif obj.type == :sphere
                        if in_sphere(x, y, z, obj.params[:center], obj.params[:radius])
                            is_solid = true
                            break
                        end
                    end
                end
                
                if is_solid
                    mask[i, j, k] = 0.0
                end
            end
        end
    end
end

end # module Geometry
