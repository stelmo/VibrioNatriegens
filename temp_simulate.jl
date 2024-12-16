using VibrioNatriegens, COBREXA, Gurobi
import AbstractFBCModels as A
import ConstraintTrees as C

label(r) = begin
    r = String(last(r))
    println(r)
    ec = join(get(model.reactions[r].annotations, "EC", [""]), ";")
    r * ": " * A.reaction_name(model, r) * (isempty(ec) ? "" : " [" * ec * "]")
end

model = VibrioNatriegens.build_model(;)
mids = A.metabolites(model)

for r in A.reactions(model)
    isnothing(A.reaction_name(model, r)) && @info r
end
for m in A.metabolites(model)
    isnothing(A.metabolite_name(model, m)) && @info m
end


println(A.n_reactions(model), ",", A.n_metabolites(model))

sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
C.pretty(C.filter_leaves(x -> abs(x) > 1e-3, sol.fluxes); format_label = r -> label(r))

