using GLM, DataFrames, COBREXA, Gurobi, AbstractFBCModels, VibrioNatriegens, ConstraintTrees, DataFramesMeta
import ConstraintTrees as C
import AbstractFBCModels as A
using CairoMakie, AlgebraOfGraphics
using CSV, JSON

function switchoff!(model, rid)
    model.reactions[rid].lower_bound = 0 
    model.reactions[rid].upper_bound = 0 
end

print_solution(sol) = C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && (startswith(string(last(ix)), "EX_") || string(last(ix)) == "BIOMASS")
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)


get_flux_prediction(sol, rid) = sol.fluxes[Symbol(rid)]

flux_tol = 1e-3
fluxes = DataFrame(CSV.File("scripts/prepared_mfa_fluxes.tab"))
@rsubset!(fluxes, sign(:LB) == sign(:UB), abs(:Flux) >= flux_tol) # must point in the same direction

model = VibrioNatriegens.build_model()
model.reactions["EX_15903"].lower_bound = 0.0 # set glucose to 0

switchoff!(model, "38215") # switch off zwf nadph isozyme
switchoff!(model, "11844") # PFL - anaerobic

# gdf = @rsubset(fluxes, :Condition == "Succinate")
# model.reactions["EX_30031"].lower_bound = -1000.0 # succinate

# gdf = @rsubset(fluxes, :Condition == "Alanine")
# model.reactions["EX_57972"].lower_bound = -1000.0 # alanine

# gdf = @rsubset(fluxes, :Condition == "Glucose")
# model.reactions["EX_15903"].lower_bound = -1000.0 # glucose

# gdf = @rsubset(fluxes, :Condition == "Ribose")
# model.reactions["EX_47013"].lower_bound = -1000.0 # ribose

gdf = @rsubset(fluxes, :Condition == "Glutamate")
model.reactions["EX_29985"].lower_bound = -14.0 # glutamate

sol = loopless_flux_balance_analysis(model, optimizer=Gurobi.Optimizer)
print_solution(sol) 

# gdf = @rsubset(fluxes, :Condition == "Acetate")
# model.reactions[""].lower_bound = -1000.0 # acetate

# gdf = @rsubset(fluxes, :Condition == "Glycerol")
# model.reactions[""].lower_bound = -1000.0 # Glycerol

cond_fluxes_dict = Dict(Symbol(r.RID) => r.Flux for r in eachrow(gdf))
mu = cond_fluxes_dict[:BIOMASS]
delete!(cond_fluxes_dict, :BIOMASS)
model.reactions["BIOMASS"].lower_bound = mu

ct = flux_balance_constraints(model)
ct *= :flux_loss^C.Constraint(; # squared relative error mean
        value = sum(
            (haskey(cond_fluxes_dict, rid) ? 1.0 : 0.0) * C.squared(get(cond_fluxes_dict, rid, 0.0) - ct.fluxes[rid].value) for
            rid in keys(ct.fluxes)
        ),
        bound = nothing,
    )

sol = COBREXA.optimized_values(ct, optimizer=Gurobi.Optimizer, objective=ct.flux_loss.value, sense=COBREXA.Minimal)


@rtransform!(gdf, :Prediction = get_flux_prediction(sol, :RID))
@rtransform!(gdf, :Difference = :Flux - :Prediction, :Relative = round(1 - :Prediction/:Flux, digits=2))

using JSON
fs = Dict(k => v for (k, v) in sol.fluxes)
open("fluxes.json", "w") do io
    JSON.print(io, fs)
end
