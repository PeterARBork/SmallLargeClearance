using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets
using Statistics
using CSV, Parameters

function loadparameters(paramsfile; varargs...)
    csvfile = CSV.read(paramsfile, NamedTuple)
    csvparams = NamedTuple{Symbol.(tuple(csvfile.symbol...))}(tuple(csvfile.value...))
    params = merge(csvparams, varargs)
    return params
end

function makeODEproblemofAD1injection(modelparams; dx = 0.01, order = 2)
    p = modelparams
    @unpack C0, σ0, x0, L, T = p
    # Parameters, variables, and derivatives
    @parameters t x
    @parameters D, v, α, Dm, Lm, Vc, Vv, Vb, S, Ql, kp
    @variables c(..) cV(..) cSAS(..) cB(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2

    membranefactor = Dm / (Lm * D)
    compartmentfactor = S * Dm / Lm

    # 1D PDE and boundary conditions
    eq = [
        Dt(c(t, x)) ~ D * Dxx(c(t, x)) - v * Dx(c(t, x)) - kp * (c(t, x) - cB(t)),
        Dt(cV(t)) ~ compartmentfactor * (c(t, 0) - α * cV(t)) / Vv - Ql * cV(t) / Vv,
        Dt(cSAS(t)) ~ compartmentfactor * (c(t, L) - α * cSAS(t)) / Vc + Ql * (cV(t) - cSAS(t)) / Vc,
        Dt(cB(t)) ~ Ql * cSAS(t) / Vb,
    ]
    bcs = [
        c(0, x) ~ C0 * exp(-(x .- x0).^2 / σ0^2),
        cSAS(0) ~ C0 * exp(-(L .- x0).^2 / σ0^2),
        cV(0) ~ C0 * exp(-(0 .- x0).^2 / σ0^2),
        cB(0) ~ 0.0,
        Dx(c(t, 0)) ~ (v/D) * c(t, 0) - membranefactor * (α * cV(t) - c(t, 0)),
        Dx(c(t, L)) ~ (v/D) * c(t, L) - membranefactor * (c(t, L) - α * cSAS(t)),
    ]

    # Space and time domains
    domains = [
        t ∈ Interval(0.0, T),
        x ∈ Interval(0.0, L),
    ]

    # PDE system
    @named pdesys = PDESystem(
        eq, bcs, domains,
        [t, x], [c(t, x), cV(t), cSAS(t), cB(t)],
        [D=>p.D, v=>p.v, α=>p.α, Dm=>p.Dm, Lm=>p.Lm, Vc=>p.Vc, Vv=>p.Vv, Vb=>p.Vb, S=>p.S, Ql=>p.Ql, kp=>p.kp],
    )

    discretization = MOLFiniteDifference([x => dx], t, approx_order=order)

    # Convert the PDE problem into an ODE problem
    prob = discretize(pdesys, discretization)

    # Solve ODE problem
    sol = solve(prob, Tsit5(), saveat = T/40)
    if checkmassconservation(sol, modelparams)
        @info "Solved AD1 problem with mass conservation"
    else
        @error "Mass not conserved in solution of AD1 problem"
    end

    return prob, sol
end

function solveAD1prob(prob, parameters; params...)
    p = merge(parameters, params)
    @parameters D, v, α, Dm, Lm, Vc, Vv, Vb, S, Ql, kp

    newp = [D=>p.D, v=>p.v, α=>p.α, Dm=>p.Dm, Lm=>p.Lm, Vc=>p.Vc, Vv=>p.Vv, Vb=>p.Vb, S=>p.S, Ql=>p.Ql, kp=>p.kp]
    newprob = remake(prob, p=newp)
    sol = solve(newprob, Tsit5(), saveat = p.T/40)

    if checkmassconservation(sol, parameters)
        return sol
    else
        return NaN
    end
end

function calcmassesovertime(sol, parameters)
    @parameters x t
    @variables c(..) cV(..) cSAS(..) cB(..)

    discrete_x = sol[x]
    solc = sol[c(t, x)]
    solcSAS = sol[cSAS(t)]
    solcV = sol[cV(t)]
    solcB = sol[cB(t)]

    @unpack S, L, Vc, Vv, Vb = parameters
    dx = sum(diff(discrete_x)) / length(discrete_x)
    mbrain = S * sum(solc * dx, dims=2)
    mSAS = Vc * solcSAS
    mV = Vv * solcV
    mB = Vb * solcB

    totalmass = mbrain .+ mSAS .+ mV .+ mB

    return Dict("total"=>totalmass, "brain"=>mbrain, "ventricle"=>mV, "SAS"=>mSAS, "blood"=>mB)
end

function checkmassconservation(sol, parameters; reltol=1e-2)
    masses = calcmassesovertime(sol, parameters)
    totalmass = masses["total"]

    coefficientofvariation = std(totalmass) / mean(totalmass)
    check = coefficientofvariation < reltol
    if check
        return true
    else
        @error "coefficient of variation is $coefficientofvariation > $reltol"
        return false
    end
end
