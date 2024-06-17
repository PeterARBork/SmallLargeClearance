using Plots, LaTeXStrings
include("../src/ADmodel.jl")
include("../src/plotfunctions.jl")

baseparameters = loadparameters("data/parametersBork2024.csv")

Dhat_Bork2024_mm2perh = 0.46
w_Bork2024_Da = 0.6 * 1e3
k1 = Dhat_Bork2024_mm2perh * w_Bork2024_Da^(1/3)
calcExpectedD_mmperh(Da) = k1 * Da^(-1/3)

vL_mm2perh = baseparameters.v * baseparameters.L
vL_mm2pers = (baseparameters.v / 60^2) * baseparameters.L

begin
    molWeights_kDa = exp.(range(log(1e-5), log(75), length=150))
    p = plot(title="flow clears large molecules", titlefontsize=12, grid=false)
    plot!(
        molWeights_kDa, x->vL_mm2perh/calcExpectedD_mmperh(x * 1e3),
        xlabel="molecular weight (kDa)",
        ylabel="bulk flow vs diffusion (Péclet)",
        label="",
    )
    plot!(3*[1.0, 1.0], [0.0, 1.0], color=:black, label="")
    plot!([0.0, 3.0], [1.0, 1.0], color=:black, label="")

    molWeights_kDa_small = exp.(range(log(1e-5), log(5), length=40))
    plot!(
        molWeights_kDa_small, x->vL_mm2perh/calcExpectedD_mmperh(x * 1e3),
        label="",
        ylabel="Péclet",
        xlabel="kDa",
        inset=(1, bbox(0., 0.23, 0.5, 0.34, :right)),
        subplot=2,
        yticks=[0, 0.5, 1.0],
        xticks=[1, 3, 5,],
    )
    plot!(3*[1.0, 1.0], [0.0, 1.0], color=:black, label="", subplot=2)
    plot!([0.0, 3.0], [1.0, 1.0], color=:black, label="", subplot=2)
end

savefig("img/weight_vs_Pe.pdf")
