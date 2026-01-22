module Visualization

using ..Common
using ..Fields
using ..Grid

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
function render_slice(
    buffers::CFDBuffers,
    grid::GridData,
    config::VizConfig,
    step::Int,
    dim_params::DimensionParams
)
    # This will be overwritten by Extension if CairoMakie is loaded
    # Or implies text output is handled here?
    # Text output doesn't need CairoMakie.
    if config.text_output
        write_slice_text("$(config.output_dir)/slice_$(step).txt", buffers, grid, config, dim_params)
    end
    # If no plotting (CairoMakie missing), specific plotting config does nothing.
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
    open(filepath, "w") do io
         println(io, "x y z u v w p")
         # Loop over slice
         # logic depending on plane...
         # ...
    end
end

end # module Visualization
