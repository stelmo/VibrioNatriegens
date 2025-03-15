using VibrioNatriegens
import AbstractFBCModels as A
using COBREXA, Gurobi
using JSON
using ConstraintTrees
import ConstraintTrees as C
using JSONFBCModels

model = VibrioNatriegens.build_model()

# model.reactions["temp"] = VibrioNatriegens.Reaction(
#     stoichiometry = Dict(
#         "29991" => -1,
#     ),
#     lower_bound = -1000.0,
#     upper_bound = 1000.0
# )
# model.reactions["temp2"] = VibrioNatriegens.Reaction(
#     stoichiometry = Dict(
#         "17319" => -1,
#     ),
#     lower_bound = -1000.0,
#     upper_bound = 1000.0
# )

model.reactions["biomass"].lower_bound = 0.6 # minimum growth rate
o2_ex = :EX_15379

ct = flux_balance_constraints(model)
sol = optimized_values(ct, optimizer=Gurobi.Optimizer, objective=ct.fluxes[o2_ex].value, sense=Maximal)
sol.fluxes[o2_ex]

# sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
# sol = parsimonious_flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

# sol = loopless_flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

open("vnat_fluxes.json", "w") do io
    JSON.print(io, sol.fluxes)
end

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && begin
            mets = [
                A.metabolite_name(model, k) for
                k in keys(A.reaction_stoichiometry(model, string(last(ix))))
            ]
            # any(in.(mets, Ref(["NADH"])))
            # any(in.(mets, Ref(["H(+)"])))
            # any(in.(mets, Ref(["Na(+)"])))
            # any(in.(mets, Ref(["O2"])))
            any(in.(mets, Ref(["H2O2"])))
            # any(in.(mets, Ref(["oxaloacetate"])))
            
            

        end
    end;
    # format_label = x -> A.reaction_name(model, string(last(x))),
)

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 100
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)
