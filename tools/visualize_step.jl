using Plots
using Statistics

include("sph_reader.jl")

function plot_bfs(step_idx)
    path = "verification/backward_step/output"
    v_file = joinpath(path, "vel_$(lpad(step_idx, 7, '0')).sph")
    p_file = joinpath(path, "prs_$(lpad(step_idx, 7, '0')).sph")
    
    if !isfile(v_file) || !isfile(p_file)
        println("File not found: $v_file or $p_file")
        return
    end

    vd = read_sph_vector(v_file)
    pd = read_sph_scalar(p_file)
    
    # Y断面 (j=23)
    j = 23
    x = [vd.x0 + (i - 0.5) * vd.dx for i in 1:vd.nx]
    z = [vd.z0 + (k - 0.5) * vd.dz for k in 1:vd.nz]
    
    u = vd.u[:, j, :]
    w = vd.w[:, j, :]
    p = pd.p[:, j, :]
    
    # 1. 速度ベクトルのみ (Quiver)
    p1 = plot(title="Velocity Vectors (Step $step_idx)", aspect_ratio=:equal, 
              xlims=(x[1], x[end]), ylims=(z[1], z[end]), xlabel="x [m]", ylabel="z [m]",
              legend=false)
    
    # 速度ベクトル (Quiver) - 間引き
    step_vec = 4
    qx = [x_val for x_val in x[1:step_vec:end], z_val in z[1:step_vec:end]]
    qz = [z_val for x_val in x[1:step_vec:end], z_val in z[1:step_vec:end]]
    qu = u[1:step_vec:end, 1:step_vec:end]
    qw = w[1:step_vec:end, 1:step_vec:end]
    
    quiver!(p1, vec(qx), vec(qz), quiver=(vec(qu), vec(qw)) .* 0.05, color=:black)
    
    # 2. 圧力等高線のみ
    p2 = contour(x, z, p', title="Pressure Contours", aspect_ratio=:equal, 
                 levels=20, color=:black, lw=0.8, xlims=(x[1], x[end]), ylims=(z[1], z[end]),
                 xlabel="x [m]", ylabel="z [m]", legend=false)
    
    # 幾何形状の追加
    for plt in (p1, p2)
        # ステップ境界
        plot!(plt, [-0.5, 0.0, 0.0], [0.3, 0.3, 0.0], color=:blue, lw=2, label="")
        # 固体領域の塗りつぶし (グレー)
        shape = Shape([-0.5, 0.0, 0.0, -0.5], [0.0, 0.0, 0.3, 0.3])
        plot!(plt, shape, color=:gray, alpha=0.3, label="")
    end
    
    final_plot = plot(p1, p2, layout=(2, 1), size=(1000, 1000))
    save_file = joinpath(path, "viz_$(step_idx)_plots.png")
    savefig(final_plot, save_file)
    println("Plot saved to $save_file")
end

# 最新のステップを確認して描画
if abspath(PROGRAM_FILE) == @__FILE__
    step_val = isempty(ARGS) ? 1000 : parse(Int, ARGS[1])
    plot_bfs(step_val)
end
