using VibrioNatriegens
import AbstractFBCModels as A
using COBREXA, Gurobi
using JSON
using ConstraintTrees
import ConstraintTrees as C
using JSONFBCModels
using JSON

model = VibrioNatriegens.build_model()

flux_balance_analysis(model,optimizer=Gurobi.Optimizer)

ex_rids = filter(startswith("EX_"), A.reactions(model))
transporters = [rid for rid in A.reactions(model) if model.reactions[rid].transporter]
for rid in [ex_rids; transporters]
    model.reactions[rid].lower_bound = -1000
    model.reactions[rid].upper_bound = 1000
end
model.reactions["ATPM"].lower_bound = 0.0

#model.reactions["biomass"].stoichiometry["597326"] = -0.5 # metabolite reaction id from escher

model.reactions["sink_reaction1"] = VibrioNatriegens.Reaction(
    stoichiometry = Dict("44841" => -1.0),
    lower_bound = -1000.0,
    upper_bound = 1000.0,
)
model.reactions["sink_reaction2"] = VibrioNatriegens.Reaction(
    stoichiometry = Dict("17071" => -1.0),
    lower_bound = -1000.0,
    upper_bound = 1000.0,
)
model.reactions["sink_reaction3"] = VibrioNatriegens.Reaction(
    stoichiometry = Dict("17001" => -1.0),
    lower_bound = -1000.0,
    upper_bound = 1000.0,
)
m = convert(JSONFBCModels.JSONFBCModel, model)
AbstractFBCModels.save(m, "vnat.json")

rid = "10543" # rhea id from blocked.csv
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

open("vnat_fluxes.json", "w") do io
    JSON.print(io, ub.fluxes)
end

C.pretty(
    C.ifilter_leaves(ub.fluxes) do ix, x
        abs(x) > 1e-6 && begin
            mets = [
                A.metabolite_name(model, k) for
                k in keys(A.reaction_stoichiometry(model, string(last(ix))))
            ]
            any(in.(mets, Ref(["pyridoxine"])))
        end
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)
