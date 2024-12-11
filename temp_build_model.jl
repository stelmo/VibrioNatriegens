using CSV, DataFrames, DataFramesMeta, XLSX, VibrioNatriegens, Gurobi
import COBREXA as X
import AbstractFBCModels as A
import ConstraintTrees as C

df = DataFrame(
    XLSX.readtable(
        joinpath("data", "curation", "curated", "base_reactions.xlsx"),
        "Sheet1",
    ),
)

rxns = VibrioNatriegens.get_kegg_reactions_in_pathway("map00010")
@rsubset!(df, :Reaction in rxns)

heteros = @rsubset(df, !ismissing(:Subunit))
homos = @rsubset(df, ismissing(:Subunit))
@select!(homos, :Reaction, :Protein, :Stoichiometry)

anno = DataFrame(CSV.File(joinpath("data", "curation", "all_annotations.csv")))
@rsubset!(anno, !ismissing(:KEGG_Reaction_Definition))
@rsubset!(anno, !ismissing(:EN_KO), !ismissing(:EN_Reaction))
@rename!(anno, :Reaction = :EN_Reaction)

@select!(
    anno,
    :Reaction,
    :Protein,
    :EN_EC,
    :EN_KO,
    :KEGG_KO,
    :HAMAP_Subunit,
    :Uniprot_Subunit,
    :EN_Symbol,
    :KEGG_Description,
    :RefSeq_Description,
    :KEGG_Reaction_Definition,
)

ghomos = groupby(homos, [:Reaction, :Protein])
gheteros = groupby(heteros, [:Reaction, :Subunit])

#! Build model

model = VibrioNatriegens.Model()

VibrioNatriegens.extend_model!(model, ghomos)
VibrioNatriegens.extend_model!(model, gheteros)
# VibrioNatriegens.curate!(model)
VibrioNatriegens.reaction_directions!(model)
VibrioNatriegens.add_exchanges!(model)
# VibrioNatriegens.add_biomass!(model)
# VibrioNatriegens.add_atpm!(model)

model
VibrioNatriegens.printmodel(model)

ct = X.parsimonious_flux_balance_constraints(model)
sol = X.optimized_values(ct; optimizer=Gurobi.Optimizer, objective=ct.objective.value)
ct.objective.bound = C.EqualTo(sol.objective)
sol = X.optimized_values(ct; optimizer=Gurobi.Optimizer, objective=ct.parsimonious_objective.value, sense=X.Minimal)

rename_func(r) = begin
    if startswith(r, "EX_")
        c = last(split(r, "_"))
        "Exchange: " * A.metabolite_name(model, string(c))
    else
        A.reaction_name(model, r)
    end
end

Dict(
    rename_func(string(k)) => v for
    (k, v) in C.filter_leaves(x -> abs(x) > 1e-3, sol.fluxes) if
    startswith(string(k), "EX_") || !startswith(string(k), "R")
)

Dict(
    rename_func(string(k)) => v for
    (k, v) in C.filter_leaves(x -> abs(x) > 1e-3, sol.fluxes) if
    startswith(string(k), "R")
)

