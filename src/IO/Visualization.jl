module Visualization

using ..Common
using ..Fields
using ..Grid
using Printf

export VizConfig, render_slice, write_slice_text

struct VizConfig
    interval::Int
    plane::Symbol         # :xy, :xz, :yz
    plane_index::Int
    variables::Vector{Symbol}  # [:velocity, :pressure]
    output_format::Symbol  # :png or :svg
    output_dir::String
    vector_enabled::Bool
    vector_skip::Int
    text_output::Bool
end

# Default implementation (stub)
# Default implementation (Generic fallback)
function render_slice(
    buffers,
    grid,
    config,
    step,
    dim_params
)
    # This will be specialized by Extension if CairoMakie is loaded.
    # The default handles text output if requested.
    if config.text_output
        write_slice_text("$(config.output_dir)/slice_$(step).txt", buffers, grid, config, dim_params)
    end
end

function write_slice_text(
    filepath::String,
    buffers::CFDBuffers,
    grid::GridData,
    config::VizConfig,
    dim_params::DimensionParams
)
    # Implement text output here as it doesn't need separate package
    mkpath(dirname(filepath))
    x = grid.x .* dim_params.L0
    y = grid.y .* dim_params.L0
    zc = grid.z_center .* dim_params.L0
    U0 = dim_params.U0

    u_dim = buffers.u .* U0
    v_dim = buffers.v .* U0
    w_dim = buffers.w .* U0
    p_dim = buffers.p .* (U0^2)

    max_x = grid.mx - 4
    max_y = grid.my - 4
    max_z = grid.mz - 4

    open(filepath, "w") do io
        println(io, "x y z u v w p")

        if config.plane == :xy
            k = clamp(config.plane_index, 1, max_z) + 2
            for j in 3:grid.my-2, i in 3:grid.mx-2
                @printf(io, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                        x[i], y[j], zc[k],
                        u_dim[i, j, k], v_dim[i, j, k], w_dim[i, j, k], p_dim[i, j, k])
            end
        elseif config.plane == :xz
            j = clamp(config.plane_index, 1, max_y) + 2
            for k in 3:grid.mz-2, i in 3:grid.mx-2
                @printf(io, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                        x[i], y[j], zc[k],
                        u_dim[i, j, k], v_dim[i, j, k], w_dim[i, j, k], p_dim[i, j, k])
            end
        elseif config.plane == :yz
            i = clamp(config.plane_index, 1, max_x) + 2
            for k in 3:grid.mz-2, j in 3:grid.my-2
                @printf(io, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                        x[i], y[j], zc[k],
                        u_dim[i, j, k], v_dim[i, j, k], w_dim[i, j, k], p_dim[i, j, k])
            end
        end
    end
end

end # module Visualization
