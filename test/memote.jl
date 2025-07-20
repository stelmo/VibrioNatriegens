
@testset "Mass balance" begin
    #=
    Check that each reaction is mass balanced.
    Note, this test also ensures that the model is stoichiometrically consistent since full mass balance implies stoichiometric consistency.
    =#
    model = deepcopy(full_model)

    rids = filter(x -> !startswith(x, "EX_") && x != "BIOMASS" && x != "BIOMASS_core", A.reactions(model))
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

@testset "Charge balance" begin
    #=
    Check that each reaction is charge balanced.
    =#
    model = deepcopy(full_model)

    rids = filter(x -> !startswith(x, "EX_") && x != "BIOMASS" && x != "BIOMASS_core", A.reactions(model))
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

@testset "ATP yield" begin
    #=
    Ensure the ATP yield per carbon source is reasonable.
    =#
    model = deepcopy(full_model)

    model.reactions["BIOMASS"].objective_coefficient = 0.0
    model.reactions["ATPM"].objective_coefficient = 1.0
    model.reactions["ATPM"].lower_bound = 0.0
    model.reactions["EX_15903"].lower_bound = 0.0 # glucose

    csources = [ # (exchange rid, atp/substrate)
        (:EX_15903, -26.0) # glucose cf. E. coli 23.5
        (:EX_57972, -11.75) # alanine  cf. E. coli 10.75
        (:EX_30031, -13.25) # succinate cf. E. coli 12.25
        (:EX_17754, -14.0) # glycerol  cf. E. coli 13.5
        (:EX_29985, -18.25) # glutamate cf. E. coli 16.5
        (:EX_47013, -20.66) # ribose  cf. E. coli -18.5
        (:EX_30089, -6.75) # acetate  cf. E. coli 6.25
        (:EX_58723, -26.0) # glucosamine  cf. E. coli 23.5
        (:EX_506227, -33.0) # n-acetyl-d-glucosamine  cf. E. coli 30
        (:EX_18391, -23.42) # gluconate  cf. E. coli 21
        (:EX_15589, -12.25) # malate  cf. E. coli 11.25
        (:EX_29806, -12.25) # fumarate  cf. E. coli 11.25
        (:EX_40886, -20.67) # L-arabinose  cf. E. coli 19.5
    ]
    ct = flux_balance_constraints(model)
    for (rid, atp) in csources
        ct.fluxes[rid].bound = C.Between(-10.0, 1000)
        sol = optimized_values(
            ct,
            optimizer = HiGHS.Optimizer,
            sense = Maximal,
            objective = ct.objective.value,
        )
        @info("ATP: $rid: $(sol.objective / sol.fluxes[rid])")
        @test isapprox(sol.objective / sol.fluxes[rid], atp, atol = 0.1)
        ct.fluxes[rid].bound = C.Between(0.0, 1000)
    end

    # anaerobic
    ct.fluxes.EX_15903.bound = C.Between(-10.0, 1000)
    ct.fluxes.EX_15379.bound = C.EqualTo(0.0)
    sol = optimized_values(
        ct,
        optimizer = HiGHS.Optimizer,
        sense = Maximal,
        objective = ct.objective.value,
    )
    @test isapprox(sol.objective / sol.fluxes.EX_15903, -3.0, atol = 1e-3) # cf. E. coli 2.5

end

@testset "ATP production" begin
    #=
    Ensure the ATP cannot be produced when the exchanges are closed.
    =#
    model = deepcopy(full_model)

    model.reactions["BIOMASS"].objective_coefficient = 0.0
    model.reactions["ATPM"].objective_coefficient = 1.0
    model.reactions["ATPM"].lower_bound = 0.0

    # test that the model cannot produce ATP with closed exchanges
    for rid in filter(x -> startswith(x, "EX_"), A.reactions(model))
        model.reactions[rid].lower_bound = 0.0
        model.reactions[rid].upper_bound = 0.0
    end
    sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    @test isapprox(sol.objective, 0.0)
end

@testset "Biomass consistency" begin
    #=
    Ensure biomass sums to 1 gram.
    =#
    model = deepcopy(full_model)

    biomass = model.reactions["BIOMASS"].stoichiometry
    btot = 0.0
    for (k, v) in biomass
        btot -=
            v * parse(Float64, first(model.metabolites[k].annotations["molarmass"])) / 1000
    end
    @test isapprox(btot, 1.0, atol = 1e-3)

    # test the core (reduced) one too
    biomass = model.reactions["BIOMASS_core"].stoichiometry
    btot = 0.0
    for (k, v) in biomass
        btot -=
            v * parse(Float64, first(model.metabolites[k].annotations["molarmass"])) / 1000
    end
    @test isapprox(btot, 1.0, atol = 1e-3)
end

@testset "Aerobic growth" begin
    #=
    Check if vibrio can be grow on carbon sources from Hoffart 2017
    =#
    model = deepcopy(full_model)
    model.reactions["EX_15903"].lower_bound = 0.0

    aerobic_substrates = [
        ("glucose", "EX_15903", 1.68),
        ("galactose", "EX_28061", 0.18),
        ("rhamnose", "EX_62346", 0.4),
        # ("maltose", "", 1.22), #! TODO: add metabolite to model
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
        # ("starch", "", 0.19), #! TODO: add metabolite to model
    ]

    for (s, e, g) in aerobic_substrates
        model.reactions[e].lower_bound = -30.0
        sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
        if s in ["maltose", "starch"]
            @test_broken sol.objective > 0.1
        else
            @test sol.objective > 0.1
        end
        model.reactions[e].lower_bound = 0.0
    end
end

@testset "Commonly secreted metabolites" begin
    #=
    Check if vibrio can secrete known products.
    =#
    model = deepcopy(full_model)
    model.reactions["BIOMASS"].lower_bound = 0.6 # minimum growth rate

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
    #=
    Check if vibrio can grow anaerobically
    =#
    model = deepcopy(full_model)
    model.reactions["EX_15379"].lower_bound = 0.0 # o2
    model.reactions["EX_15903"].lower_bound = -37.0 # glc (measured)
    sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
    @test_broken abs(1 - sol.objective / 0.57) <= 0.2 # growth rel to measured value by me (0.92 in other papers)
end

@testset "Known growth supporting substrates" begin
    #=
    Check if vibrio can be grow on carbon sources from Coppens 2023
    =#
    df = DataFrame(CSV.File(joinpath("experiments", "coppens_2023_biolog.csv")))
    dropmissing!(df)
    @select!(df, :Substrate, :Experiment, :Chebi)
    @transform!(df, :Exchange = "EX_" .* string.(:Chebi))
    @subset!(df, :Experiment)

    model = deepcopy(full_model)
    model.reactions["EX_15903"].lower_bound = 0.0

    for (s, ex) in zip(df.Substrate, df.Exchange)
        model.reactions[ex].lower_bound = -30.0
        sol = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
        if s == "Formic Acid"
            @test_broken sol.objective > 0.1
        else
            @test sol.objective > 0.1
        end
        model.reactions[ex].lower_bound = 0.0
    end
end

@testset "Growth rate tests" begin
    #=
    Check that vibrio's predicted growth rates are close to experimentally measured values when constraining the substrates and products.
    =#
    measurements = Dict(
        "alanine" => Dict(
            "BIOMASS" => ("BIOMASS", 0.91),
            "acetate" => ("EX_30089", 4.5),
            "alanine" => ("EX_57972", -30.03),
        ),
        "ribose" =>
            Dict("BIOMASS" => ("BIOMASS", 0.87), "ribose" => ("EX_47013", -16.1)),
        "glucose" => Dict(
            "BIOMASS" => ("BIOMASS", 1.7),
            "glucose" => ("EX_15903", -25.0),
            "acetate" => ("EX_30089", 14.1),
            "succinate" => ("EX_30031", 0.18),
        ),
        "glutamate" =>
            Dict("BIOMASS" => ("BIOMASS", 0.58), "glutamate" => ("EX_29985", -15.3)),
        "glycerol" =>
            Dict("BIOMASS" => ("BIOMASS", 0.62), "glycerol" => ("EX_17754", -16.5)),
        "succinate" => Dict(
            "BIOMASS" => ("BIOMASS", 1.07),
            "succinate" => ("EX_30031", -28.3),
            "fumarate" => ("EX_29806", 0.48),
        ),
    )

    model = deepcopy(full_model)
    model.reactions["EX_15903"].lower_bound = 0.0

    for k in keys(measurements)
        ct = flux_balance_constraints(model)
        for (kk, vv) in measurements[k]
            kk == "BIOMASS" && continue
            ct.fluxes[Symbol(vv[1])].bound = C.Between(vv[2], vv[2])
        end
        sol = optimized_values(
            ct,
            optimizer = HiGHS.Optimizer,
            sense = Maximal,
            objective = ct.objective.value,
        )
        mu = measurements[k]["BIOMASS"][2]
        @info "$k: $mu vs $(sol.objective)"
        @test abs(1 - sol.objective / mu) <= 0.5
    end
end
