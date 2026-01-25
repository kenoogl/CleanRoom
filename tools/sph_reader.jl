"""
SPH Reader Utility for V-Isio format
"""

function read_sph_header(io::IO)
    # Rec 1: svType, dType
    read(io, Int32)
    svType = read(io, Int32)
    dType = read(io, Int32)
    read(io, Int32)
    
    # Rec 2: nx, ny, nz
    read(io, Int32)
    nx = read(io, Int32)
    ny = read(io, Int32)
    nz = read(io, Int32)
    read(io, Int32)
    
    # Rec 3: time, step, origin, pitch
    read(io, Int32)
    time = read(io, Float32)
    step = read(io, Int32)
    x0 = read(io, Float32)
    y0 = read(io, Float32)
    z0 = read(io, Float32)
    dx = read(io, Float32)
    dy = read(io, Float32)
    dz = read(io, Float32)
    read(io, Int32)
    
    return (nx=nx, ny=ny, nz=nz, x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, time=time, step=step, svType=svType)
end

function read_sph_vector(filepath::String)
    open(filepath, "r") do io
        h = read_sph_header(io)
        n = h.nx * h.ny * h.nz
        bytes_per = n * sizeof(Float32)

        # Rec 4: Data (either interleaved or block-separated)
        rec_bytes = read(io, Int32)
        if rec_bytes == 3 * bytes_per
            # Interleaved u,v,w per cell
            data = Vector{Float32}(undef, 3 * n)
            read!(io, data)
            read(io, Int32)

            u = zeros(Float32, h.nx, h.ny, h.nz)
            v = zeros(Float32, h.nx, h.ny, h.nz)
            w = zeros(Float32, h.nx, h.ny, h.nz)
            idx = 1
            for k in 1:h.nz, j in 1:h.ny, i in 1:h.nx
                u[i, j, k] = data[idx]
                v[i, j, k] = data[idx+1]
                w[i, j, k] = data[idx+2]
                idx += 3
            end
            return (u=u, v=v, w=w, nx=h.nx, ny=h.ny, nz=h.nz, x0=h.x0, y0=h.y0, z0=h.z0, dx=h.dx, dy=h.dy, dz=h.dz, time=h.time, step=h.step)
        elseif rec_bytes == bytes_per
            # Block-separated u, v, w (writer in this repo)
            data_u = Vector{Float32}(undef, n)
            read!(io, data_u)
            read(io, Int32)

            rec2 = read(io, Int32)
            if rec2 != bytes_per
                error("Invalid SPH record size for v: expected $(bytes_per) got $(rec2)")
            end
            data_v = Vector{Float32}(undef, n)
            read!(io, data_v)
            read(io, Int32)

            rec3 = read(io, Int32)
            if rec3 != bytes_per
                error("Invalid SPH record size for w: expected $(bytes_per) got $(rec3)")
            end
            data_w = Vector{Float32}(undef, n)
            read!(io, data_w)
            read(io, Int32)

            u = reshape(data_u, h.nx, h.ny, h.nz)
            v = reshape(data_v, h.nx, h.ny, h.nz)
            w = reshape(data_w, h.nx, h.ny, h.nz)
            return (u=u, v=v, w=w, nx=h.nx, ny=h.ny, nz=h.nz, x0=h.x0, y0=h.y0, z0=h.z0, dx=h.dx, dy=h.dy, dz=h.dz, time=h.time, step=h.step)
        else
            error("Unsupported SPH record size: $(rec_bytes)")
        end
    end
end

function read_sph_scalar(filepath::String)
    open(filepath, "r") do io
        h = read_sph_header(io)
        
        # Rec 4: Data
        read(io, Int32)
        data = Vector{Float32}(undef, h.nx * h.ny * h.nz)
        read!(io, data)
        read(io, Int32)
        p = reshape(data, h.nx, h.ny, h.nz)
        return (p=p, nx=h.nx, ny=h.ny, nz=h.nz, x0=h.x0, y0=h.y0, z0=h.z0, dx=h.dx, dy=h.dy, dz=h.dz, time=h.time, step=h.step)
    end
end
