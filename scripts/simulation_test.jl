using VibrioNatriegens
import AbstractFBCModels as A
using COBREXA, Gurobi
using JSON
using ConstraintTrees
import ConstraintTrees as C
using JSONFBCModels

model = VibrioNatriegens.build_model()

model.reactions["EX_15903"].lower_bound = 0.0
model.reactions["EX_47013"].lower_bound = -15.0

sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

sol = parsimonious_flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
sol = loopless_flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

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
            any(in.(mets, Ref(["NADH"])))
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




