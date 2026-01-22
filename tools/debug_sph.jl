function inspect_sph_header(filepath)
    open(filepath, "r") do io
        # Rec 1
        r1 = read(io, Int32)
        sv = read(io, Int32)
        dt = read(io, Int32)
        read(io, Int32)
        println("Rec1: size=$r1 svType=$sv dType=$dt")
        
        # Rec 2
        r2 = read(io, Int32)
        nx = read(io, Int32)
        ny = read(io, Int32)
        nz = read(io, Int32)
        read(io, Int32)
        println("Rec2: size=$r2 nx=$nx ny=$ny nz=$nz")
        
        # Rec 3
        r3 = read(io, Int32)
        time = read(io, Float32)
        step = read(io, Int32)
        x0 = read(io, Float32)
        y0 = read(io, Float32)
        z0 = read(io, Float32)
        dx = read(io, Float32)
        dy = read(io, Float32)
        dz = read(io, Float32)
        read(io, Int32)
        println("Rec3: size=$r3")
        println("  time=$time step=$step")
        println("  x0=$x0 y0=$y0 z0=$z0")
        println("  dx=$dx dy=$dy dz=$dz")
    end
end

inspect_sph_header(ARGS[1])
