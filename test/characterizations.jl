
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
