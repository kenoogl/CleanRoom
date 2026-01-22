module SPHWriter

using ..Common
using ..Grid

export generate_sph_filename, write_sph_vector, write_sph_scalar

function generate_sph_filename(prefix::String, step::Int)::String
    return "$(prefix)_$(lpad(step, 7, '0')).sph"
end

function write_record(io::IO, data...)
    # Calculate size bytes
    total_bytes = 0
    for d in data
        total_bytes += sizeof(d)
    end
    
    write(io, Int32(total_bytes))
    for d in data
        write(io, d)
    end
    write(io, Int32(total_bytes))
end

function write_array_record(io::IO, arr::Array{Float32})
    total_bytes = sizeof(arr)
    write(io, Int32(total_bytes))
    write(io, arr)
    write(io, Int32(total_bytes))
end

function write_sph_header(
    io::IO,
    svType::Int32,
    dType::Int32,
    grid::GridData,
    step::Int,
    time::Float64,
    is_dim::Bool,
    dim_params::DimensionParams
)
    # Rec 1: svType, dType
    write_record(io, svType, dType)
    
    # Rec 2: IM, JM, KM
    IM = Int32(grid.mx-4) # Internal size? V-Isio usually visualizes internal.
    JM = Int32(grid.my-4)
    KM = Int32(grid.mz-4)
    # But we calculate on full grid including ghosts?
    # Usually visualization handles physical domain.
    # Outputting full grid with ghosts is easier, but user might want physical domain.
    # Task 13.4 (Checkpoint) says "ゴーストセル込み".
    # SPH is for Viz. Usually we crop ghosts.
    # "GridData" stores mx, my, mz (inc ghosts).
    # x indices 3 to mx-2 are physical. size = mx-4.
    # I'll output Internal Only.
    
    nx, ny, nz = grid.mx-4, grid.my-4, grid.mz-4
    write_record(io, Int32(nx), Int32(ny), Int32(nz))
    
    # Rec 3: Time, Step, Org, Pitch
    # Time (Float32)
    t_val = Float32(time * (is_dim ? 1.0 : dim_params.T0))
    stp = Int32(step)
    
    # Origin and Pitch
    # Origin should be at the left edge (face) of the first cell, not cell center
    # Cell center = origin + (i - 0.5) * dx
    # grid.x[3] is cell center of first internal cell
    # Left edge = grid.x[3] - dx/2 = grid.origin[1] (for internal domain)
    x0 = Float32((grid.x[3] - grid.dx/2) * dim_params.L0)
    y0 = Float32((grid.y[3] - grid.dy/2) * dim_params.L0)
    z0 = Float32((grid.z_center[3] - grid.dz[3]/2) * dim_params.L0)
    dx = Float32(grid.dx * dim_params.L0)
    dy = Float32(grid.dy * dim_params.L0)
    dz = Float32(grid.dz[3] * dim_params.L0) # Using first cell dz
    
    write_record(io, t_val, stp, x0, y0, z0, dx, dy, dz)
end

function write_sph_vector(
    filepath::String,
    u::Array{Float64, 3},
    v::Array{Float64, 3},
    w::Array{Float64, 3},
    grid::GridData,
    step::Int,
    time::Float64,
    is_dimensional::Bool,
    dim_params::DimensionParams
)
    # Crop to internal domain
    # u[3:mx-2, 3:my-2, 3:mz-2]
    mx, my, mz = grid.mx, grid.my, grid.mz
    nx, ny, nz = mx-4, my-4, mz-4
    
    U_scale = is_dimensional ? 1.0 : dim_params.U0
    
    u_out = Float32.(u[3:mx-2, 3:my-2, 3:mz-2] .* U_scale)
    v_out = Float32.(v[3:mx-2, 3:my-2, 3:mz-2] .* U_scale)
    w_out = Float32.(w[3:mx-2, 3:my-2, 3:mz-2] .* U_scale)
    
    open(filepath, "w") do io
        write_sph_header(io, Int32(2), Int32(1), grid, step, time, is_dimensional, dim_params)
        
        # Interleave u, v, w? Or u block, v block, w block?
        # V-Isio Vector standard: u, v, w interleaved (u1, v1, w1, u2, ...) ?
        # Or Block?
        # "Fortran unformatted" usually implies one huge record array.
        # Often (u, v, w) separate records or one record (u,v,w).
        # Common convention: u, then v, then w.
        # Checking similar codes (H2 reference): Block format (u array, v array, w array).
        # But if svType=2 (Vector), it might expect vector per point.
        # I'll write u, v, w as separate records or one combined?
        # Safest "One Record" containing u,v,w concatenated.
        # Or (u, v, w) interleaved.
        # Let's write separate records for safety or interleaved?
        # Interleaved is common for "Vector Field".
        # But Arrays in Julia are column major.
        # I'll create a single buffer of size 3*N.
        
        data = Vector{Float32}(undef, 3 * nx * ny * nz)
        idx = 1
        # Loop order x, y, z (Julia standard)
        for k in 1:nz, j in 1:ny, i in 1:nx
             data[idx]   = u_out[i,j,k]
             data[idx+1] = v_out[i,j,k]
             data[idx+2] = w_out[i,j,k]
             idx += 3
        end
        write_array_record(io, data)
    end
end

function write_sph_scalar(
    filepath::String,
    p::Array{Float64, 3},
    grid::GridData,
    step::Int,
    time::Float64,
    is_dimensional::Bool,
    dim_params::DimensionParams
)
    mx, my, mz = grid.mx, grid.my, grid.mz
    nx, ny, nz = mx-4, my-4, mz-4
    
    P_scale = is_dimensional ? 1.0 : dim_params.U0^2 # Approx
    
    p_out = Float32.(p[3:mx-2, 3:my-2, 3:mz-2] .* P_scale)
    
    open(filepath, "w") do io
        write_sph_header(io, Int32(1), Int32(1), grid, step, time, is_dimensional, dim_params)
        write_array_record(io, vec(p_out))
    end
end

end # module SPHWriter
