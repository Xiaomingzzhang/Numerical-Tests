using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, GLMakie
using GLMakie.GeometryBasics
f1(x, p, t) = SA[x[2], p[1]*x[1]+p[3]*sin(2pi * t)]
f2(x, p, t) = SA[x[2], -p[2]*x[1]+p[3]*sin(2pi * t)]
hyper1(x, p, t) = x[1] - (p[4]+p[5]*cos(2pi * t))
hyper2(x, p, t) = x[1] + (p[4]+p[5]*cos(2pi * t))
dom1(x, p, t) = -p[4]-p[5]*cos(2pi * t) < x[1]
dom2(x, p, t) = x[1] < -(p[4]+p[5]*cos(2pi * t))
impact_rule(x, p, t) = SA[x[1], -p[6]*x[2]-(1+p[6])*p[5]*2*pi*sin(2pi * t)]
id(x,p,t) = x
vectorfield = PiecewiseImpactV((f1, f2), (dom1, dom2), (hyper1, hyper2), (impact_rule, id), [1])
setup = setmap(vectorfield, (0.0, 1.0), Tsit5(), abstol=1e-8, reltol=1e-8, cross_time=0.0001)
para = [2, 5, 0.5, 2, 0.1, 0.98]
df1(x, p, t) = SA[0 1; p[1] 0]
initialguess = SA[0.0, 0.0]
saddle = findsaddle(f1, df1, (0.0, 1.0), initialguess, para, abstol=1e-10)
iscontact(setup, saddle, para)
prob = NSOneDManifoldProblem(setup, para)
segment = gen_segment(saddle)
manifold = growmanifold(prob, segment, 7)
manifold
function manifold_plot(data)
    set_theme!(theme_latexfonts(), fontsize=22)
    fig = Figure(size=(800, 640))
    axes = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", xlabelsize=26, ylabelsize=26, ylabelrotation=false, xticks=range(-2.1,2.1,length=5))
    for k in eachindex(data)
        for j in eachindex(data[k])
            points=Point2f.(data[k][j].u)
            lines!(axes,points)
        end
    end
    fig
end
function manifold_plot_no_ticks(data)
    set_theme!(theme_latexfonts(), fontsize=22)
    fig = Figure(size=(800, 640))
    axes = Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", xlabelsize=26, ylabelsize=26, ylabelrotation=false)
    for k in eachindex(data)
        for j in eachindex(data[k])
            points=Point2f.(data[k][j].u)
            lines!(axes,points)
        end
    end
    fig
end
fig1 = manifold_plot(manifold.data)
fig2 = manifold_plot_no_ticks(manifold.data)
