module Checkpoint

using ..Common
using ..Fields
using ..Grid

export generate_checkpoint_filename, write_checkpoint, read_checkpoint

function generate_checkpoint_filename(step::Int)::String
    return "checkpoint_$(lpad(step, 7, '0')).bin"
end

function write_checkpoint(
    filepath::String,
    buffers::CFDBuffers,
    grid::GridData,
    step::Int,
    time::Float64,
    dim_params::DimensionParams
)
    Sys.islittleendian() || error("Checkpoint requires little-endian system.")
    open(filepath, "w") do io
        # Header
        # Spec: internal counts, Int32 header
        write(io, Int32(grid.mx - 4))
        write(io, Int32(grid.my - 4))
        write(io, Int32(grid.mz - 4))
        write(io, Int32(step))
        write(io, Float64(time))
        write(io, Int32(0)) # is_dimensional (always 0)
        write(io, Float64(dim_params.L0))
        write(io, Float64(dim_params.U0))
        
        # Arrays (Raw write of memory)
        # u, v, w, p
        write(io, buffers.u)
        write(io, buffers.v)
        write(io, buffers.w)
        write(io, buffers.p)
    end
end

function read_checkpoint(
    filepath::String,
    buffers::CFDBuffers,
    dim_params::DimensionParams
)::Tuple{Int, Float64}
    if !isfile(filepath)
        error("Checkpoint file not found: $filepath")
    end
    Sys.islittleendian() || error("Checkpoint requires little-endian system.")
    
    step = 0
    time = 0.0
    
    open(filepath, "r") do io
        nx = read(io, Int32)
        ny = read(io, Int32)
        nz = read(io, Int32)
        if nx != size(buffers.u, 1) - 4 || ny != size(buffers.u, 2) - 4 || nz != size(buffers.u, 3) - 4
            error("Checkpoint grid size mismatch: File($nx,$ny,$nz) vs Current($(size(buffers.u, 1)-4),$(size(buffers.u, 2)-4),$(size(buffers.u, 3)-4))")
        end
        
        step = Int(read(io, Int32))
        time = read(io, Float64)
        is_dim = read(io, Int32)
        L0 = read(io, Float64)
        U0 = read(io, Float64)
        
        if abs(L0 - dim_params.L0) > 1e-6 || abs(U0 - dim_params.U0) > 1e-6
            println("Warning: Checkpoint dim params mismatch (L0: $L0 vs $(dim_params.L0), U0: $U0 vs $(dim_params.U0))")
        end
        if is_dim != 0
            error("Checkpoint file is not nondimensional (is_dimensional=$(is_dim)).")
        end
        
        read!(io, buffers.u)
        read!(io, buffers.v)
        read!(io, buffers.w)
        read!(io, buffers.p)
    end
    
    return (step, time)
end

end # module Checkpoint
