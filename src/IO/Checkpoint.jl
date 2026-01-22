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
    is_dimensional::Bool,
    dim_params::DimensionParams
)
    open(filepath, "w") do io
        # Header
        # Nx, Ny, Nz (Internal counts used for grid gen, or array sizes?)
        # Better save Array Sizes mx, my, mz.
        # But restart needs to rebuild GridConfig/Fields?
        # Usually restart loads params JSON first, rebuilds grid/fields, then OVERWRITES with checkpoint data.
        # So we verify sizes match.
        
        write(io, grid.mx)
        write(io, grid.my)
        write(io, grid.mz)
        write(io, step)
        write(io, time)
        write(io, is_dimensional)
        write(io, dim_params.L0)
        write(io, dim_params.U0)
        
        # Arrays (Raw write of memory)
        # u, v, w, p
        write(io, buffers.u)
        write(io, buffers.v)
        write(io, buffers.w)
        write(io, buffers.p)
        
        # Save avg if count > 0?
        # For full restartability of averaging.
        write(io, buffers.avg_count[])
        write(io, buffers.u_avg)
        write(io, buffers.v_avg)
        write(io, buffers.w_avg)
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
    
    step = 0
    time = 0.0
    
    open(filepath, "r") do io
        mx = read(io, Int)
        my = read(io, Int)
        mz = read(io, Int)
        
        if mx != size(buffers.u, 1) || my != size(buffers.u, 2) || mz != size(buffers.u, 3)
            error("Checkpoint grid size mismatch: File($mx,$my,$mz) vs Current($(size(buffers.u)))")
        end
        
        step = read(io, Int)
        time = read(io, Float64)
        is_dim = read(io, Bool)
        L0 = read(io, Float64)
        U0 = read(io, Float64)
        
        if abs(L0 - dim_params.L0) > 1e-6 || abs(U0 - dim_params.U0) > 1e-6
            println("Warning: Checkpoint dim params mismatch (L0: $L0 vs $(dim_params.L0), U0: $U0 vs $(dim_params.U0))")
        end
        
        read!(io, buffers.u)
        read!(io, buffers.v)
        read!(io, buffers.w)
        read!(io, buffers.p)
        
        cnt = read(io, Int)
        buffers.avg_count[] = cnt
        read!(io, buffers.u_avg)
        read!(io, buffers.v_avg)
        read!(io, buffers.w_avg)

        # Dimensionality check
        if is_dim
            # If stored as dimensional, must nondimensionalize
            # u_nd = u / U0
            # p_nd = p / (rho*U0^2) ? rho=1 in code?
            # Checkpoint spec 13.5 says "有次元データの無次元化処理".
            # Assuming rho=1 for now.
            # u, v, w /= U0
            # p /= U0^2 ? (Pressure scaling)
            # x, y, z are grid, not stored in checkpoint.
            
            # Note: Current implementation writes raw buffers.
            # `write_checkpoint` takes `is_dimensional`.
            # If caller says it's dimensional, it implies buffers contain dimensional data?
            # But buffers usually store non-dimensional data during run.
            # So is_dimensional usually False in checkpoint unless specified.
            # If restart source is dimensional, convert:
            
            U_scale = U0
            P_scale = U0^2
            
            @. buffers.u /= U_scale
            @. buffers.v /= U_scale
            @. buffers.w /= U_scale
            @. buffers.p /= P_scale
            # Also avgs
            @. buffers.u_avg /= U_scale
            @. buffers.v_avg /= U_scale
            @. buffers.w_avg /= U_scale
        end
    end
    
    return (step, time)
end

end # module Checkpoint
