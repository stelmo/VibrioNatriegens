using GLM, DataFrames, COBREXA, Gurobi, AbstractFBCModels, VibrioNatriegens, ConstraintTrees
import ConstraintTrees as C
import AbstractFBCModels as A
using CairoMakie, AlgebraOfGraphics

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

# split the difference between 13C MFA and actual measurements
measurements = Dict(
    "alanine" => Dict("acetate" => ("EX_30089", 6.3), "alanine" => ("EX_57972", -35.03), "BIOMASS" => ("BIOMASS_core", 0.95)),
    "ribose" => Dict("ribose" => ("EX_47013", -16.1), "BIOMASS" => ("BIOMASS_core", 0.871)),
    "glucose" => Dict(
        "glucose" => ("EX_15903", -21.34),
        "acetate" => ("EX_30089", 15.5),
        "succinate" => ("EX_30031", 0.19),
        "BIOMASS" => ("BIOMASS_core", 1.696)
    ),
    "glutamate" => Dict("glutamate" => ("EX_29985", -14.9), "BIOMASS" => ("BIOMASS_core", 0.576)),
    "glycerol" => Dict("glycerol" => ("EX_17754", -13.5), "BIOMASS" => ("BIOMASS_core", 0.634), "glutamate" => ("EX_29985", 0.27)),
    "succinate" => Dict("succinate" => ("EX_30031", -23.64), "fumarate" => ("EX_29806", 0.52), "BIOMASS" => ("BIOMASS_core", 1.074)),
    "acetate" => Dict("acetate" => ("EX_30089", -26.6), "BIOMASS" => ("BIOMASS_core", 0.424)),
    "salt" => Dict(
        "glucose" => ("EX_15903", -10.766),
        "acetate" => ("EX_30089", 8.83),
        "succinate" => ("EX_30031", 0.22),
        "BIOMASS" => ("BIOMASS_core", 0.595)
    ),
    "iptg" => Dict("glucose" => ("EX_15903", -27.16), "acetate" => ("EX_30089", 19.52), "BIOMASS" => ("BIOMASS_core", 1.897)),
)

mus = Dict(k => v["BIOMASS"][2] for (k, v) in measurements)

model = VibrioNatriegens.build_model()
model.reactions["EX_15903"].lower_bound = 0.0 # set glucose to 0
model.reactions["BIOMASS"].objective_coefficient = 0.0
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


draw(
    data(df) * mapping(:mu, :maxatp) * (visual(Scatter) + linear())
)

# using CSV
# CSV.write("atp_fit.csv", df)
