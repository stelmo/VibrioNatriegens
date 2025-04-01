using GLM, DataFrames, COBREXA, Gurobi, AbstractFBCModels, VibrioNatriegens, ConstraintTrees
import ConstraintTrees as C
import AbstractFBCModels as A
using CairoMakie

#=
Predict appropriate GAM & NGAM values using the techniques in https://www.nature.com/articles/nbt.3956

Importantly, before running the code below, set the GAM and NGAM to zero in the model build scripts!
This can be found in `src/biomass.jl`.
=#

print_solution(sol) = C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)

measurements = Dict(
    "alanine" => Dict("acetate" => ("EX_30089", 4.5), "alanine" => ("EX_57972", -40.03), "biomass" => ("biomass", 0.916)),
    "ribose" => Dict("ribose" => ("EX_47013", -16.1), "biomass" => ("biomass", 0.871)),
    "glucose" => Dict(
        "glucose" => ("EX_15903", -25.0),
        "acetate" => ("EX_30089", 14.1),
        "succinate" => ("EX_30031", 0.18),
        "biomass" => ("biomass", 1.696)
    ),
    "glutamate" => Dict("glutamate" => ("EX_29985", -15.3), "biomass" => ("biomass", 0.576)),
    "glycerol" => Dict("glycerol" => ("EX_17754", -16.5), "biomass" => ("biomass", 0.634)),
    "succinate" => Dict("succinate" => ("EX_30031", -28.3), "fumarate" => ("EX_29806", 0.48), "biomass" => ("biomass", 1.074)),
    "acetate" => Dict("acetate" => ("EX_30089", -26.6), "biomass" => ("biomass", 0.424)),
    "salt" => Dict(
        "glucose" => ("EX_15903", -10.766),
        "acetate" => ("EX_30089", 8.83),
        "succinate" => ("EX_30031", 0.22),
        "biomass" => ("biomass", 0.595)
    ),
    "iptg" => Dict("glucose" => ("EX_15903", -27.16), "acetate" => ("EX_30089", 19.52), "biomass" => ("biomass", 1.897)),
)

mus = Dict(
    "iptg" => 1.897,
    "salt" => 0.595,
    "alanine" => 0.916,
    "acetate" => 0.424,
    "ribose" => 0.871,
    "glucose" => 1.696,
    "glutamate" => 0.576,
    "glycerol" => 0.634,
    "succinate" => 1.074,
)

model = VibrioNatriegens.build_model()
model.reactions["EX_15903"].lower_bound = 0.0 # set glucose to 0
model.reactions["biomass"].objective_coefficient = 0.0
model.reactions["ATPM"].objective_coefficient = 1.0
model.reactions["ATPM"].lower_bound = 0.0

df = DataFrame(maxatp = Float64[], mu = Float64[], id=String[])
for k in keys(measurements)
    ct = flux_balance_constraints(model)
    for (kk, vv) in measurements[k]
        ct.fluxes[Symbol(vv[1])].bound = C.EqualTo(vv[2])
    end
    sol = optimized_values(
        ct,
        optimizer = Gurobi.Optimizer,
        sense = Maximal,
        objective = ct.objective.value,
    )
    # print_solution(sol)
    push!(df, (sol.objective, mus[k], k))
end

ols = lm(@formula(maxatp ~ mu), df)
df
a, b = coef(ols)


fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Growth rate [1/h]", ylabel = "Max ATP flux [mmol/gDW/h]")
ablines!(ax, [a], [b])
scatter!(ax, df.mu, df.maxatp)
for r in eachrow(df)
    text!(ax, r.mu, r.maxatp, text=r.id)
end
fig
