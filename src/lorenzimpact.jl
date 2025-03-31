using InvariantManifolds, LinearAlgebra, StaticArrays, OrdinaryDiffEq, GLMakie,DataInterpolations
function lorenz(x, p, t)
    σ, ρ, β = p
    v = SA[σ*(x[2]-x[1]),
        ρ*x[1]-x[2]-x[1]*x[3],
        x[1]*x[2]-β*x[3]
    ]
    v / sqrt(0.1 + norm(v)^2)
end
hyper(x,p,t) = x[3]-p[4]
rule(x,p,t) = SA[x[1], -x[2], x[3]]
vectorfield =BilliardV(lorenz, (hyper,),(rule,))
setup = setmap(vectorfield, (0.0, -1.0), Tsit5(), abstol=1e-8)
para = [10.0, 28.0, 8/3, 10.0]
function eigenv(p)
    σ, ρ, β = p
    [SA[0.0, 0.0, 1.0], SA[-(-1 + σ + sqrt(1 - 2 * σ + 4 * ρ * σ + σ^2))/(2*ρ), 1, 0]]
end
saddle = Saddle(SA[0, 0, 0.0], eigenv(para), [1.0, 1.0])
disk = gen_disk(saddle, r=1.0, d=0.1)
prob = NSVTwoDManifoldProblem(setup, para, amax=0.5, d=0.5, ϵ=0.2, dsmin=1e-3)
manifold = growmanifold(prob, disk, 120, interp=LinearInterpolation)
manifold
function manifold_plot(annulus)
    fig = Figure()
    axes = LScene(fig[1, 1], show_axis=false, scenekw=(backgroundcolor=:white, clear=true))
    second(x) = x[2]
    for i in eachindex(annulus)
        for j in eachindex(annulus[i])
            points = annulus[i][j].u
            lines!(axes, first.(points), second.(points), last.(points), fxaa=true)
        end
    end
    fig
end
fig = manifold_plot(manifold.data)

