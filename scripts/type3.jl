
using VibrioNatriegens
import AbstractFBCModels as A
using COBREXA, Gurobi
using JSON
using ConstraintTrees
import ConstraintTrees as C
using JSONFBCModels


model = VibrioNatriegens.build_model()

sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)

for rid in filter(x -> startswith(x, "EX_"), A.reactions(model))
    model.reactions[rid].lower_bound = 0.0
    model.reactions[rid].upper_bound = 0.0
end

model.reactions["ATPM"].lower_bound = 0.0
m = flux_balance_constraints(model)

fva = Dict()
for rid in A.reactions(model)
    ub = optimized_values(
        m,
        optimizer = Gurobi.Optimizer,
        objective = m.fluxes[Symbol(rid)].value,
        sense = COBREXA.Maximal,
    ).objective
    lb = optimized_values(
        m,
        optimizer = Gurobi.Optimizer,
        objective = m.fluxes[Symbol(rid)].value,
        sense = COBREXA.Minimal,
    ).objective
    fva[rid] = (lb, ub)
end
fva
