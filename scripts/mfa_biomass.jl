using GLM, DataFrames, COBREXA, Gurobi, AbstractFBCModels, VibrioNatriegens, ConstraintTrees, DataFramesMeta
import ConstraintTrees as C
import AbstractFBCModels as A
using CairoMakie, AlgebraOfGraphics
using CSV

#=
Predict appropriate GAM & NGAM values using the techniques in https://www.nature.com/articles/nbt.3956

Importantly, before running the code below, set the GAM to zero in the model build scripts!
This can be found in `src/biomass.jl`.
=#

extendbound(lb, ub, fact) = begin
   _lb = lb - abs(lb) * fact
   _ub = ub + abs(ub) * fact
   
   if lb >= 0 && _lb < 0
        _lb = 0
   end
   
   if ub <= 0 && _ub > 0 
        _ub = 0
   end

   (_lb, _ub)
end

flux_tol = 1e-3
fluxes = DataFrame(CSV.File("scripts/prepared_mfa_fluxes.tab"))
@rsubset!(fluxes, sign(:LB) == sign(:UB), abs(:Flux) >= flux_tol) # must point in the same direction
@rsubset!(fluxes, :Condition != "Glycerol")

model = VibrioNatriegens.build_model()
model.reactions["EX_15903"].lower_bound = 0.0 # set glucose to 0
model.reactions["BIOMASS"].objective_coefficient = 0.0
model.reactions["ATPM"].objective_coefficient = 1.0
model.reactions["ATPM"].lower_bound = 0.0

# switch off zwf nadph isozyme
model.reactions["38215"].lower_bound = 0 
model.reactions["38215"].upper_bound = 0 

# PFL - anaerobic
model.reactions["11844"].lower_bound = 0 
model.reactions["11844"].upper_bound = 0 

df = DataFrame(maxatp = Float64[], mu = Float64[], id = String[])
for gdf in groupby(fluxes, :Condition) 

    cond = first(gdf.Condition)
    @info cond

    cond_fluxes_dict = Dict(Symbol(r.RID) => (r.LB, r.UB) for r in eachrow(gdf))
    (mulb, muub) = cond_fluxes_dict[:BIOMASS]
    delete!(cond_fluxes_dict, :BIOMASS)

    ct = flux_balance_constraints(model)
    ct.fluxes.BIOMASS.bound = C.Between(mulb, muub)
    
    for (k, (lb, ub)) in cond_fluxes_dict
        if startswith(string(k), "EX_")
            ct.fluxes[k].bound = C.Between(lb, ub) # fix 
        else
           ct.fluxes[k].bound = C.Between(extendbound(lb, ub, 0.5)...)
        end
    end

    sol = optimized_values(
        ct,
        optimizer = Gurobi.Optimizer,
        sense = Maximal,
        objective = ct.objective.value,
    )

    # print_solution(sol)
    if isnothing(sol)
        @warn cond
    else
        push!(df, (sol.objective, sol.fluxes.BIOMASS, cond))
    end
end

ols = lm(@formula(maxatp ~ mu), df)
df
a, b = coef(ols)

draw(
    data(df) * mapping(:mu, :maxatp) * (visual(Scatter) + mapping(text = :id => verbatim)*visual(Annotation)) + 
    data(df) * mapping(:mu, :maxatp) * (visual(Scatter) + linear()),
    axis = (limits = ((0, 2.0), (0, 250)),),
)

# using CSV
# CSV.write("atp_fit.csv", df)

# m = deepcopy(model)
# k = "glutamate"
# for (kk, vv) in measurements[k]
#     m.reactions[vv[1]].upper_bound = vv[2]
#     m.reactions[vv[1]].lower_bound = vv[2]
# end

# sol = parsimonious_flux_balance_analysis(m, optimizer=Gurobi.Optimizer)

# sol = loopless_flux_balance_analysis(m; optimizer = Gurobi.Optimizer)

# using JSON
# fs = Dict(k => v for (k, v) in sol.fluxes)
# open("fluxes.json", "w") do io
#     JSON.print(io, fs)
# end
