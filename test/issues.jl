@testset "Known blocked reactions" begin
    #=
    These are known blocked reactions that could not be resolved - gaps in knowledge
    =#
    blocked_reactions = [
        "10412"
        "14265"
        "15805"
        "16861"
        "17029"
        "17417"
        "17769"
        "20712"
        "22060"
        "22472"
        "22748"
        "24360"
        "27485"
        "30815"
        "30935"
        "31115"
        "32527"
        "33867"
        "34095"
        "34111"
        "42288"
        "42700"
        "77643"
        "nfn"
    ]

    model = VibrioNatriegens.build_model()
    exchange_rids = filter(startswith("EX_"), A.reactions(model))
    for rid in exchange_rids
        model.reactions[rid].lower_bound = -1000
        model.reactions[rid].upper_bound = 1000
    end
    model.reactions["ATPM"].lower_bound = 0.0

    vs = screen(blocked_reactions) do rid
        ct = flux_balance_constraints(model)
        lb = optimized_values(
            ct,
            objective = ct.fluxes[rid].value,
            sense = Minimal,
            optimizer = HiGHS.Optimizer,
        )
        ub = optimized_values(
            ct,
            objective = ct.fluxes[rid].value,
            sense = Maximal,
            optimizer = HiGHS.Optimizer,
        )
        max(
            isnothing(lb) ? 0 : abs(lb.fluxes[rid]),
            isnothing(ub) ? 0 : abs(ub.fluxes[rid]),
        )
    end

    for v in vs
        @test_broken v != 0
    end
end
