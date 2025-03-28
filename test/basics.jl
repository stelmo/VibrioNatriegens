
@testset "Mass balanced" begin
    # this test also ensures that the model is stoichiometrically consistent IF all the reactions mass balance.
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

@testset "Sensible ATP production" begin

    model = VibrioNatriegens.build_model()

    model.reactions["biomass"].objective_coefficient = 0.0
    model.reactions["ATPM"].objective_coefficient = 1.0
    model.reactions["ATPM"].lower_bound = 0.0

    sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    atp_per_glucose = abs(sol.fluxes.ATPM / sol.fluxes.EX_15903)
    @test isapprox(atp_per_glucose, 32.25, atol = 1e-3) # e coli -23.5

    for rid in filter(x -> startswith(x, "EX_"), A.reactions(model))
        model.reactions[rid].lower_bound = 0.0
        model.reactions[rid].upper_bound = 0.0
    end

    sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    @test isapprox(sol.objective, 0.0)

end

@testset "Biomass consistency" begin

    model = VibrioNatriegens.build_model()

    biomass = model.reactions["biomass"].stoichiometry
    btot = 0.0
    for (k, v) in biomass
        btot -=
            v * parse(Float64, first(model.metabolites[k].annotations["molarmass"])) / 1000
    end

    @test isapprox(btot, 1.0, atol = 1e-3)
end

@testset "Aerobic growth" begin
    # based on Hoffart 2017
    model = VibrioNatriegens.build_model()
    model.reactions["EX_15903"].lower_bound = 0.0

    aerobic_substrates = [
        ("glucose", "EX_15903", 1.68),
        ("galactose", "EX_28061", 0.18),
        ("rhamnose", "EX_62346", 0.4),
        ("maltose", "EX_15903", 1.22), #! this is glucose, maltose is not a metabolite yet TODO
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
        ("starch", "EX_15903", 0.19), #! this is glucose, starch is not a metabolite yet TODO
    ]

    for (s, e, g) in aerobic_substrates
        model.reactions[e].lower_bound = -20.0
        sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
        # println(s, ": ", g, " vs ", sol.objective)
        # push!(gs, (s, sol.objective > 0.1))
        if s in ["maltose", "starch"]
            @test_broken sol.objective > 0.1
        else
            @test sol.objective > 0.1
        end
        model.reactions[e].lower_bound = 0.0
    end
end

@testset "Commonly secreted metabolites" begin

    model = VibrioNatriegens.build_model()
    model.reactions["biomass"].lower_bound = 0.6 # minimum growth rate

    secreted_products = [
        :EX_30089 # acetate
        :EX_30031 # succinate
        :EX_15740 # formate
        :EX_29806 # fumarate
        :EX_16004 # lactate
        :EX_16236 # ethanol
        :EX_57762 # valine
        :EX_29985 # glutamate
        :EX_57416 # alanine
        :EX_15361 # pyruvate
        :EX_16526 # CO2
    ]

    ct = flux_balance_constraints(model)
    for ex_rid in secreted_products
        sol = optimized_values(
            ct,
            optimizer = HiGHS.Optimizer,
            objective = ct.fluxes[ex_rid].value,
            sense = Maximal,
        )
        @test sol.fluxes[ex_rid] > 5.0 # minimum flux
    end
end

@testset "Anaerobic growth" begin

    model = VibrioNatriegens.build_model()
    model.reactions["EX_15379"].lower_bound = 0.0
    sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    @test sol.objective > 0.1
    # @test abs(1 - sol.objective / 0.92) <= 0.3 # growth rel to measured value
end

@testset "Known growth supporting substrates" begin

    df = DataFrame(CSV.File(joinpath("experiments", "coppens_2023_biolog.csv")))
    dropmissing!(df)
    @select!(df, :Substrate, :Experiment, :Chebi)
    @transform!(df, :Exchange = "EX_" .* string.(:Chebi))
    @subset!(df, :Experiment)

    model = VibrioNatriegens.build_model()
    model.reactions["EX_15903"].lower_bound = 0.0

    for (s, ex) in zip(df.Substrate, df.Exchange)
        model.reactions[ex].lower_bound = -10.0
        sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
        res = isnothing(sol) ? false : (sol.objective > 0.1)
        if s == "Formic Acid"
            @test_broken res
        else
            res || @warn(s)
            @test res
        end
        model.reactions[ex].lower_bound = 0.0
    end
end
