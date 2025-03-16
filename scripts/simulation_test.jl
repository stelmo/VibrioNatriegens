using VibrioNatriegens
import AbstractFBCModels as A
using COBREXA, Gurobi
using JSON
using ConstraintTrees
import ConstraintTrees as C
using JSONFBCModels

model = VibrioNatriegens.build_model()

sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
sol = parsimonious_flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

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
