using Plots, LaTeXStrings
include("../src/ADmodel.jl")
include("../src/plotfunctions.jl")

baseparameters = loadparameters("data/parametersBork2024.csv")
parameters = merge(
    baseparameters,
    (
        T=10,
        C0=10,
        σ0=0.4,
        x0=0.5,
        D=baseparameters.D,
        Dm=baseparameters.D,
        #v=0,
        L=3,
    )
)
#Miaodata = loaddata("Miaodata.csv")

prob, baselinesol = makeODEproblemofAD1injection(parameters);
baselinesol = solveAD1prob(prob, parameters)
lowDfac = 5
baselinesol_lowD = solveAD1prob(prob, parameters; D=parameters.D / lowDfac, Dm=parameters.Dm / lowDfac)
Qlfac = 10
highQlsol = solveAD1prob(prob, parameters; Ql=parameters.Ql * Qlfac)
highQlsol_lowD = solveAD1prob(prob, parameters; Ql=parameters.Ql * Qlfac, D=parameters.D / lowDfac, Dm=parameters.D / lowDfac)

begin
    timepoint = 3
    shuntx = 0.9 * parameters.L
    p = plot(layout=(2, 2), size=(700, 400),)

    plot!(
        title="brain concentration at $(timepoint)h", titlefontsize=12,
        xlabel="", ylabel=L"c(%$timepoint, x)",
        subplot=1, legend=:bottom,
    )
    plotsoloverspace!(baselinesol; color=2, #ylims=(0, 1.2),
     timepoint, label=L"KX $(D_m=D)$", subplot=1)
    plotsoloverspace!(highQlsol; color=:black, timepoint, label=L"fast shunt $(Q_l\times %$Qlfac)$", subplot=1)

    plot!(title="", xlabel=L"$x$ (mm)", ylabel=L"c(%$timepoint, x)", subplot=3, titlefontsize=12)
    plotsoloverspace!(baselinesol, timepoint=0.0; color=:gray, linestyle=:dash, label="init.", subplot=3)
    plotsoloverspace!(baselinesol_lowD; color=2, timepoint, label="KX, slow D", subplot=3)
    plotsoloverspace!(highQlsol_lowD; color=:black, timepoint, label="fast shunt, slow D", subplot=3)

    plot!(title="fiber model", titlefontsize=12, xlabel="", ylabel=L"c(t, %$shuntx)", subplot=2)
    plotmodelfiber!(baselinesol, shuntx; color=2, label="", subplot=2)
    plotmodelfiber!(highQlsol, shuntx, ylims=(0, 1.2); color=:black, label="", subplot=2)

    plotmodelfiber!(baselinesol_lowD, shuntx; color=2, label="",ylabel=L"c(t, %$shuntx)", subplot=4)
    plotmodelfiber!(highQlsol_lowD, shuntx, ylims=(0, 1.2); color=:black, label="", subplot=4, xlabel=L"time $t$ (h)",)
end
savefig("img/kx_vs_high_Ql.pdf")

highk_p_sol = baselinesol_lowD = solveAD1prob(prob, parameters; k_p=1.0)

@parameters t x
@variables c(..) cV(..) cSAS(..) cB(..)

sol = highk_p_sol
discrete_x = sol[x]
discrete_t = sol[t]
solc = sol[c(t, x)]
solcV = sol[cV(t,)]
solcSAS = sol[cSAS(t,)]
solcB = sol[cB(t,)]

baselinesolc = baselinesol[c(t, x)]
baselinesolcV = baselinesol[cV(t,)]
baselinesolcSAS = baselinesol[cSAS(t,)]
baselinesolcB = baselinesol[cB(t,)]

discrete_x = sol[x]
discrete_t = sol[t]

anim = @animate for i in 1:length(discrete_t)
    p1 = bar([1., ], [solcV[i], ], ylim=[0, 10], xticks=([1, ], ["ventricle"]), label="", ylabel="c", title="ventricle")
    bar!([1., ], [baselinesolcV[i], ], ylim=[0, 10], xticks=([1, ], ["ventricle"]), label="", ylabel="c", title="ventricle")

    p2 = plot(discrete_x, solc[i, :], label=""; legend=false, title="brain", xlabel="mm", ylim=[0,10])
    plot!(discrete_x, baselinesolc[i, :], label=""; legend=false, title="brain", xlabel="mm", ylim=[0,10])

    p3 = bar([1, 2, ], [solcSAS[i], solcB[i], ], xticks=([1, 2], ["SAS", "blood"]), ylim=[0, 10], label="", title="SAS, blood")
    bar!([1, 2, ], [baselinesolcSAS[i], baselinesolcB[i], ], xticks=([1, 2], ["SAS", "blood"]), ylim=[0, 10], label="", title="SAS, blood")

    plot(p1, p2, p3, layout=(1, 3))
end
gif(anim, "high_kp.gif",fps=10)
