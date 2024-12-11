
function reaction_directions!(model; excludes = String[], threshold = 3.0)

    rxns = filter(x -> x ∉ excludes, collect(keys(model.reactions)))
    filter!(x -> !_is_cached("directionality", x), rxns)

    if !isempty(rxns)
        eq = eQuilibrator.Equilibrator() # slow, only do if necessary
        for rid in rxns
            println(rid)
            rxnstoich = model.reactions[rid].stoichiometry
            substrates = [string(abs(Int(v))) * " " * k for (k, v) in rxnstoich if v < 0]
            products = [string(abs(Int(v))) * " " * k for (k, v) in rxnstoich if v > 0]
            rxn_string = kegg(join(substrates, " + ") * " = " * join(products, " + "))
            _d = ln_reversibility_index(eq, rxn_string; skip_unbalanced = true)
            _pdg = physiological_dg_prime(eq, rxn_string; skip_unbalanced = true)
            pdg = isnothing(_pdg) ? nothing : Measurements.value(ustrip(u"kJ/mol", _pdg))
            d = isnothing(_d) ? nothing : Measurements.value(ustrip(u"kJ/mol", _d))
            _cache("directionality", rid, (d, pdg))
        end
    end

    for rid in keys(model.reactions)

        (dg, dgg) = VibrioNatriegens._get_cache("directionality", rid)
        model.reactions[rid].dg = dgg

        if isnothing(dg) || (abs(dg) <= threshold)
            model.reactions[rid].lower_bound = -1000
            model.reactions[rid].upper_bound = 1000
        elseif dg < -threshold # forward
            model.reactions[rid].lower_bound = 0
            model.reactions[rid].upper_bound = 1000
        elseif dg > threshold # reverse
            model.reactions[rid].lower_bound = -1000
            model.reactions[rid].upper_bound = 0
        else
            @warn("Directions confusion!")
        end
    end

end
