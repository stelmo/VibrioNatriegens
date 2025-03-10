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

model.reactions["biomass"].stoichiometry["597326"] = -0.5 # pyridoxal
model.reactions["biomass"].stoichiometry["57287"] = -0.5 # coa

#model.reactions["sink_reaction1"] = VibrioNatriegens.Reaction( # need this
#   stoichiometry = Dict("57287" => -1.0),
#   lower_bound = 0.0,
#   upper_bound = 1000.0,
#)

rid = "10227" # rhea id from blocked.csv 
#rid = "sink_reaction1"
ct = flux_balance_constraints(model)
sol1 = optimized_values(
    ct,
    objective = ct.fluxes[rid].value,
    sense = Minimal,
    optimizer = Gurobi.Optimizer,
)
sol2 = optimized_values(
    ct,
    objective = ct.fluxes[rid].value,
    sense = Maximal,
    optimizer = Gurobi.Optimizer,
)
lb = sol1.fluxes[rid]
ub = sol2.fluxes[rid]

ct *= :parsimony^C.Constraint(
    sum(C.squared(x.value) for x in values(ct.fluxes)),
    nothing
)
ubct = deepcopy(ct)
ubct *= :ub^C.Constraint(
    ubct.fluxes[Symbol(rid)].value,
    C.EqualTo(ub),
)
sol2 = optimized_values(
    ubct,
    objective = ct.parsimony.value,
    sense = Minimal,
    optimizer = Gurobi.Optimizer,
)

open("vnat_fluxes.json", "w") do io
    JSON.print(io, sol2.fluxes)
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
