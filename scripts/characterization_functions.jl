
function blocked_reactions()
    # open all exchanges find reactions that cannot carry flux under any circumstances
    @everywhere begin
        model = VibrioNatriegens.build_model()
        exchange_rids = filter(startswith("EX_"), A.reactions(model))
        for rid in exchange_rids
            model.reactions[rid].lower_bound = -1000
            model.reactions[rid].upper_bound = 1000
        end
        model.reactions["ATPM"].lower_bound = 0.0
    end

    non_exchange_rids = setdiff(A.reactions(model), exchange_rids)
    vs = screen(non_exchange_rids, workers = workers()) do rid
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

    length(non_exchange_rids[vs .== 0])
end

function deadend_metabolites()
    #=
    A dead-end metabolite (DEM) is defined as a metabolite that is produced by
    the known metabolic reactions of an organism and has no reactions consuming
    it, or that is consumed by the metabolic reactions of an organism and has no
    known reactions producing it, and in both cases has no identified
    transporter
    =#

end
