using VibrioNatriegens
import AbstractFBCModels as A
using COBREXA, Gurobi
using JSON
using ConstraintTrees
import ConstraintTrees as C
using JSONFBCModels, JSON

block_rxn(model, rid) = begin
    model.reactions[rid].lower_bound = 0.0
    model.reactions[rid].upper_bound = 0.0
end
block_rxns(model, rids) = begin
    for rid in rids
        block_rxn(model, string(rid))
    end
end

model = VibrioNatriegens.build_model()
model.reactions["EX_15379"].lower_bound = 0.0 #
model.reactions["EX_15379"].upper_bound = 0.0 #
model.reactions["EX_15903"].lower_bound = -35.0 #

# ct.fluxes.EX_15903.bound = C.Between(-10.0, 1000)
# ct.fluxes.EX_15379.bound = C.EqualTo(0.0)    

# # these are involved in loops...
# block_rxns(model, [
#     23148
#     28046
#     28046
#     30699
#     30699
#     19125
#     22852
# ])

sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
sol = parsimonious_flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

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
            # any(in.(mets, Ref(["NADH", "NAD(+)"])))
            # any(in.(mets, Ref(["H(+)"])))
            # any(in.(mets, Ref(["Na(+)"])))
            # any(in.(mets, Ref(["phosphoenolpyruvate"])))
            any(in.(mets, Ref(["pyruvate"])))
            # any(in.(mets, Ref(["O2"])))
            # any(in.(mets, Ref(["ATP"])))
            # any(in.(mets, Ref(["oxaloacetate"])))
        end
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 100
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)

#####################
using JSONFBCModels, AbstractFBCModels, COBREXA
model = convert(A.CanonicalModel.Model, COBREXA.load_model("iML1515.json"))
model.reactions["EX_glc__D_e"].lower_bound = -25.0
model.reactions["BIOMASS_Ec_iML1515_core_75p37M"].objective_coefficient = 0.0
model.reactions["ATPM"].objective_coefficient = 1.0

sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
sol = parsimonious_flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

sol.fluxes.ATPM / sol.fluxes.EX_glc__D_e # 23.5

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && begin
            mets = collect(keys(A.reaction_stoichiometry(model, string(last(ix)))))
            any(in.(mets, Ref(["nadh_c"])))
            # any(in.(mets, Ref(["atp_c"])))
        end
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)

############################
