
function plotsoloverspace!(sol; timepoint=NaN, plotvarargs...)
    @parameters t x
    @variables c(..)
    discrete_x = sol[x]
    discrete_t = sol[t]

    solc = sol[c(t, x)]

    if isnan(timepoint)
        numTs = 5
        showTindices = round.(Int, exp.(range(log(1), log(length(discrete_t)), numTs)))
        for i in showTindices
            plot!(discrete_x, solc[i, :], label = L"t=%$(discrete_t[i])"; plotvarargs...)
        end
    else
        Tindex = argmax(timepoint .== discrete_t)
        timepoint = discrete_t[Tindex]
        plot!(discrete_x, solc[Tindex, :], label=L"t=%$timepoint"; plotvarargs...)
    end
end

function plotsolovertime!(sol; plotvarargs...)
    @parameters t x
    @variables c(..) cV(..) cSAS(..) cB(..)

    discrete_x = sol[x]
    discrete_t = sol[t]

    solc = sol[c(t, x)]

    numlocs = 5
    showlocindices = round.(Int, range(1, length(discrete_x), numlocs))

    for i in showlocindices
        plot!(discrete_t, solc[:, i], label = L"x=%$(discrete_x[i])"; plotvarargs...)
    end

    solcSAS = sol[cSAS(t)]
    solcV = sol[cV(t)]
    solcB = sol[cB(t)]

    plot!(discrete_t, solcSAS, label="SAS", color=:darkblue; plotvarargs...)
    plot!(discrete_t, solcV, label="vent.", color=:lightblue; plotvarargs...)
    plot!(discrete_t, solcB, label="blood", color=:red; plotvarargs...)
end

function plotmodelfiber!(sol, fiberlocation;  plotvarargs...)
    @parameters t x
    @variables c(..)
    discrete_x = sol[x]
    discrete_t = sol[t]
    solc = sol[c(t, x)]
    fiberindex = argmax(fiberlocation .== discrete_x)
    xfiber = discrete_x[fiberindex]
    maxy = maximum(solc[:, :])
    plot!(discrete_t, solc[:, fiberindex], label=L"c(t, %$xfiber)", ylims=(0, maxy); plotvarargs...)
end

function plotconcsandmassses(sol, parameters)
    masses = calcmassesovertime(sol, parameters)
    discrete_t = sol[t]

    plot(layout=(1, 3))
    plotsoloverspace!(sol, subplot=1)
    plotsolovertime!(sol, subplot=2, legend=:topright)
    plot!(discrete_t, masses["total"], label="total", ylabel="mass", subplot=3)
    plot!(discrete_t, masses["brain"], label="brain", subplot=3)
    plot!(discrete_t, masses["ventricle"], label="vent", subplot=3)
    plot!(discrete_t, masses["SAS"], label="SAS", subplot=3)
    plot!(discrete_t, masses["blood"], label="blood", subplot=3)
end
