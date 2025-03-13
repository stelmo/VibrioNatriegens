
@testset "Mass balanced" begin

    model = VibrioNatriegens.build_model()

    rids = filter(x -> !startswith(x, "EX_") && x != "biomass", A.reactions(model))
    unbal_rids = String[]
    for rid in rids
        s = A.reaction_stoichiometry(model, rid)
        m = Dict()
        for (k, v) in s
            for (kk, vv) in A.metabolite_formula(model, k)
                m[kk] = get(m, kk, 0) + vv * v
            end
        end
        m
        all(values(m) .== 0) || push!(unbal_rids, rid)
    end

    @test isempty(unbal_rids)
end

@testset "Charge balanced" begin

    model = VibrioNatriegens.build_model()

    rids = filter(x -> !startswith(x, "EX_") && x != "biomass", A.reactions(model))
    unbal_rids = String[]
    for rid in rids
        s = A.reaction_stoichiometry(model, rid)
        m = 0
        for (k, v) in s
            m += v * A.metabolite_charge(model, k)
        end
        m
        m == 0 || push!(unbal_rids, rid)
    end

    @test isempty(unbal_rids)
end

@testset "ATP production" begin

    model = VibrioNatriegens.build_model()

    model.reactions["biomass"].objective_coefficient = 0.0
    model.reactions["ATPM"].objective_coefficient = 1.0
    model.reactions["ATPM"].lower_bound = 0.0


    sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
    llsol = loopless_flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
    @test isapprox(sol.objective, llsol.objective)

    for rid in filter(x -> startswith(x, "EX_"), A.reactions(model))
        model.reactions[rid].lower_bound = 0.0
        model.reactions[rid].upper_bound = 0.0
    end

    sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
    @test isapprox(sol.objective, 0.0)

end

@testset "Biomass consistency" begin

    model = VibrioNatriegens.build_model()

    biomass = model.reactions["biomass"].stoichiometry
    btot = 0.0
    for (k, v) in biomass
        btot = v * model.metabolites[k].molarmass
    end

    @test isapprox(btot, 1.0, atol = 1e-3)
end

@testset "Deadend metabolites" begin

    # open all exchanges, transporters, and find metabolites that can only be produced or consumed
    @everywhere begin
        model = VibrioNatriegens.build_model()
        ex_rids = filter(startswith("EX_"), A.reactions(model))
        transporters =
            [rid for rid in A.reactions(model) if model.reactions[rid].transporter]
        for rid in [ex_rids; transporters]
            model.reactions[rid].lower_bound = -1000
            model.reactions[rid].upper_bound = 1000
        end
        model.reactions["ATPM"].lower_bound = 0.0
        model.reactions["temp_reaction"] = VibrioNatriegens.Reaction(
            stoichiometry = Dict(),
            lower_bound = -1000.0,
            upper_bound = 1000.0,
        )
    end

    mids = filter(!occursin("_"), A.metabolites(model))

    vs = screen(mids, workers = workers()) do mid
        model.reactions["temp_reaction"].stoichiometry = Dict(mid => -1.0)
        ct = flux_balance_constraints(model)

        lb = optimized_values(
            ct,
            objective = ct.fluxes.temp_reaction.value,
            sense = Minimal,
            optimizer = Gurobi.Optimizer,
        )
        ub = optimized_values(
            ct,
            objective = ct.fluxes.temp_reaction.value,
            sense = Maximal,
            optimizer = Gurobi.Optimizer,
        )
        (
            isnothing(lb) ? 0 : abs(lb.fluxes.temp_reaction),
            isnothing(ub) ? 0 : abs(ub.fluxes.temp_reaction),
        )
    end

end

@testset "Blocked reactions" begin
    # open all exchanges, transporters, and find reactions that cannot carry flux under any circumstances
    @everywhere begin
        model = VibrioNatriegens.build_model()
        ex_rids = filter(startswith("EX_"), A.reactions(model))
        transporters =
            [rid for rid in A.reactions(model) if model.reactions[rid].transporter]
        for rid in [ex_rids; transporters]
            model.reactions[rid].lower_bound = -1000
            model.reactions[rid].upper_bound = 1000
        end
        model.reactions["ATPM"].lower_bound = 0.0
    end

    met_rids = setdiff(A.reactions(model), [ex_rids; transporters])
    vs = screen(met_rids, workers = workers()) do rid
        ct = flux_balance_constraints(model)
        lb = optimized_values(
            ct,
            objective = ct.fluxes[rid].value,
            sense = Minimal,
            optimizer = Gurobi.Optimizer,
        )
        ub = optimized_values(
            ct,
            objective = ct.fluxes[rid].value,
            sense = Maximal,
            optimizer = Gurobi.Optimizer,
        )
        max(
            isnothing(lb) ? 0 : abs(lb.fluxes[rid]),
            isnothing(ub) ? 0 : abs(ub.fluxes[rid]),
        )
    end

    met_rids[vs.==0]
    using CSV
    CSV.write("blocked.csv", DataFrame(Blocked = met_rids[vs.==0]))
end

@testset "Aerobic growth rates" begin
    # based on Hoffart 2017
    model = VibrioNatriegens.build_model()
    model.reactions["EX_15903"].lower_bound = 0.0

    aerobic_substrates = [
        ("glucose", "EX_15903", 1.68),
        ("galactose", "EX_28061", 0.18),
        ("rhamnose", "EX_62346", 0.4),
        # ("maltose", "", 1.22),
        ("arabinose", "EX_46994", 0.83),
        ("glycerol", "EX_17754", 0.86),
        ("fructose", "EX_28645", 1.51),
        ("sucrose", "EX_17992", 1.79),
        ("nacetylglucosamine", "EX_506227", 1.74),
        ("acetate", "EX_30089", 0.45),
        ("malate", "EX_15589", 0.85),
        ("fumarate", "EX_29806", 0.99),
        ("succinate", "EX_30031", 1.0),
        ("gluconate", "EX_18391", 1.51),
        # ("starch", "", 0.19),
    ]
    gs = []
    for (s, e, g) in aerobic_substrates
        model.reactions[e].lower_bound = -10.0
        sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
        push!(gs, (s, sol.objective > 0.1))
        model.reactions[e].lower_bound = 0.0
    end
    all(last.(gs))
end

@testset "Anaerobic growth" begin

    model = VibrioNatriegens.build_model()
    model.reactions["EX_15379"].lower_bound = 0.0
    sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

    @test abs(1 - sol.objective / 0.92) <= 0.2 # growth
    @test abs(sol.fluxes["EX_30089"]) > 0.05 # acetate
    @test abs(sol.fluxes["EX_30031"]) > 0.05 # succinate
    @test abs(sol.fluxes["EX_15740"]) > 0.05 # formate
    @test abs(sol.fluxes["EX_16004"]) > 0.05 # lactate
    @test abs(sol.fluxes["EX_16236"]) > 0.05 # ethanol
    @test abs(sol.fluxes["EX_57762"]) > 0.05 # valine
    @test abs(sol.fluxes["EX_29985"]) > 0.05 # glutamate
    @test abs(sol.fluxes["EX_57416"]) > 0.05 # alanine
end

@testset "Known growth supporting substrates" begin

    df = DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "experiments", "coppens_2023_biolog.csv"),
        ),
    )
    dropmissing!(df)
    @select!(df, :Substrate, :Experiment, :Chebi)
    @transform!(df, :Exchange = "EX_" .* last.(split.(:Chebi, ":")))
    @subset!(df, :Experiment)

    model = VibrioNatriegens.build_model()
    model.reactions["EX_15903"].lower_bound = 0.0

    res = Bool[]
    for (s, ex) in zip(df.Substrate, df.Exchange)
        println(s)
        model.reactions[ex].lower_bound = -10.0
        sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
        push!(res, isnothing(sol) ? false : (sol.objective > 0.1))
        model.reactions[ex].lower_bound = 0.0
    end
    @test all(res)

end
