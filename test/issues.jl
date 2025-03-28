@testset "Known issues" begin

    blocked_reactions = [
        "10412"
        "11196"
        "11492"
        "11528"
        "13325"
        "13541"
        "14265"
        "15613"
        "15689"
        "15805"
        "15837"
        "16861"
        "16881"
        "16913"
        "17029"
        "17105"
        "17417"
        "17753"
        "17769"
        "17801"
        "17805"
        "18985"
        "20217"
        "20517"
        "20549"
        "20712"
        "21152"
        "21708"
        "21896"
        "22024"
        "22060"
        "22472"
        "22712"
        "22748"
        "23256"
        "23560"
        "23596"
        "24360"
        "24927"
        "25796"
        "27345"
        "27485"
        "28342"
        "30119"
        "30131"
        "30367"
        "30815"
        "30935"
        "30939"
        "31107"
        "31115"
        "32459"
        "32527"
        "33707"
        "33867"
        "34095"
        "34099"
        "34111"
        "34115"
        "35095"
        "40663"
        "42288"
        "42700"
        "45332"
        "54828"
        "70223"
        "72067"
        "72075"
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
