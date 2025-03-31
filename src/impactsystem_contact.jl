using LinearAlgebra, StaticArrays, OrdinaryDiffEq, GLMakie, InvariantManifolds
f(x, p, t) = SA[x[2], sin(x[1])+p[1]*x[2]+p[2]*sin(p[3]*t)]
hyper1(x, p, t) = x[1] - p[4]
hyper2(x, p, t) = x[1] + p[4]
rule(x, p, t) = SA[x[1], -p[5]*x[2]]
vectorfield = BilliardV(f, (hyper1, hyper2), (rule,rule))
para = [0.1, 0.1, 1.64, pi/4 , 0.9]
setup = setmap(vectorfield, (0.0, 2*pi/para[3]), Tsit5(), abstol=1e-8, reltol=1e-8)
saddle = findsaddle(setup, SA[0.33105711548504757, -0.2795863847274121], para)
solver = ns_solver(vectorfield, (0.0, 2*pi/para[3]), Tsit5(), 2, Float64; abstol=1e-8, reltol=1e-8, dtmax=0.01)
sol = solver(saddle.saddle,para).sol
prob = NSOneDManifoldProblem(setup, para,d=0.001,ϵ=1e-8, dsmin=1e-6)
seg = gen_segment(saddle)
manifold = growmanifold(prob, seg, 10)
manifold
function manifold_plot(data)
    set_theme!(theme_latexfonts(), fontsize=22)
    fig = Figure(size=(800, 640))
    axes = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", xlabelsize=26, ylabelsize=26, ylabelrotation=false, xticks=([-pi / 4, -pi / 8, 0, pi / 8, pi / 4], [L"-\pi/4", L"-\pi/8", L"0", L"\pi/8", L"\pi/4"]))
    for k in eachindex(data)
        for j in eachindex(data[k])
            points = data[k][j].u
            lines!(axes, first.(points), last.(points))
        end
    end
    fig
end

function period_solution_plot(data)
    xs = data[1,:]
    ys = data[2,:]
    set_theme!(theme_latexfonts(), fontsize=22)
    fig = Figure(size=(800, 640))
    axes = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", xlabelsize=26, ylabelsize=26, ylabelrotation=false, xticks=([pi / 16, pi / 8, 3pi / 16, pi / 4], [L"\pi/16", L"\pi/8", L"3\pi/16", L"\pi/4"]))
    lines!(axes, xs, ys, color=:black)
    fig
end

function manifold_plot_no_ticks(data)
    set_theme!(theme_latexfonts(), fontsize=22)
    fig = Figure(size=(800, 640))
    axes = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", xlabelsize=26, ylabelsize=26, ylabelrotation=false)
    lines!(axes, [pi / 4, pi / 4], [-2, 2], color=:black)
    for k in eachindex(data)
        for j in eachindex(data[k])
            points = data[k][j].u
            lines!(axes, first.(points), last.(points))
        end
    end
    fig
end
fig1 = period_solution_plot(sol)
fig2 = manifold_plot(manifold.data)
fig3 = manifold_plot_no_ticks(manifold.data)