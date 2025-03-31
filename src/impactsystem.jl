using InvariantManifolds, StaticArrays, OrdinaryDiffEq, GLMakie, DataInterpolations
f(x, p, t) = SA[x[2], sin(x[1])+p[1]*x[2]+p[2]*sin(p[3] * t)]
hyper1(x, p, t) = x[1] + p[4]
hyper2(x, p, t) = x[1] - p[4]
rule(x, p, t) = SA[x[1], -p[5]*x[2]]
vectorfield = BilliardV(f, (hyper1, hyper2), (rule, rule))
para = [0.1, 0.1, 2pi, pi / 4, 0.9]
setup = setmap(vectorfield, (0.0, 2pi / para[3]), Tsit5(),
        abstol=1e-8, reltol=1e-8)
df(x, p, t) = SA[0 1; cos(x[1]) p[1]]
initialguess = SA[0.0, 0.0]
saddle = findsaddle(f, df, (0.0, 2pi / para[3]), initialguess, para, abstol=1e-10)
prob = NSOneDManifoldProblem(setup, para, dsmin=1e-7)
segment = gen_segment(saddle)
manifold = growmanifold(prob, segment, 13)
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

function manifold_plot_no_ticks(data)
    set_theme!(theme_latexfonts(), fontsize=22)
    fig = Figure(size=(800, 640))
    axes = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", xlabelsize=26, ylabelsize=26, ylabelrotation=false)
    lines!(axes, [-pi / 4, -pi / 4], [-2, 2], color=:black)
    for k in eachindex(data)
        for j in eachindex(data[k])
            points = data[k][j].u
            lines!(axes, first.(points), last.(points))
        end
    end
    fig
end
fig1 = manifold_plot(manifold.data)
fig2 = manifold_plot_no_ticks(manifold.data)