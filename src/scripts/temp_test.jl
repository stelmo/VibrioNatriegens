using VibrioNatriegens, COBREXA, Gurobi
import AbstractFBCModels as A
import ConstraintTrees as C

label(r) = begin
    r = String(last(r))
    ec = join(get(model.reactions[r].annotations, "EC", [""]), ";")
    r * ": " * A.reaction_name(model, r) * (isempty(ec) ? "" : " [" * ec * "]")
end

model = VibrioNatriegens.build_model(;)

# sol = flux_balance_analysis(model, optimizer = Gurobi.Optimizer)
# C.pretty(C.filter_leaves(x -> abs(x) > 1e-3, sol.fluxes); format_label = r -> label(r))

ignorelist = [
    "R00959",
    "R00710",
    "R00711",
    "R00754",
    "R00405",
    "R01541",
    "R08572",
    "R01741",
    "R01542",
    "R00264",
    "R02167",
    "R00271",
    "R01213",
    "R00742",
    "R00921",
    "R01353",
    "R03045",
    "R01360",
    "R00994",
    "R03896",
    "R03898",
    "R01395", # nucleotides
    "R01288", # cysteine
    "R01777", # cysteine
    "R01931", # cysteine
]
rxns = filter(startswith("R"), A.reactions(model))
filter!(x -> x ∉ ignorelist, rxns)

out = filter(
    x -> typeof(x) != Bool,
    screen(rxns) do rxn
        c = flux_balance_constraints(model)
        c.objective.value = c.fluxes[Symbol(rxn)].value
        v = max(
            abs(
                optimized_values(
                    c,
                    objective = c.objective.value,
                    optimizer = Gurobi.Optimizer,
                    output = c.objective,
                    sense = Minimal,
                ),
            ),
            abs(
                optimized_values(
                    c,
                    objective = c.objective.value,
                    optimizer = Gurobi.Optimizer,
                    output = c.objective,
                    sense = Maximal,
                ),
            ),
        )
        println(rxn)
        v > 1.0 || label((:xxx, Symbol(rxn)))
    end,
)

println(A.n_reactions(model), ",", A.n_metabolites(model))


using VibrioNatriegens
VibrioNatriegens.parse_rhea_reactions(joinpath("data", "rhea", "rhea-reactions.txt"))
VibrioNatriegens.parse_chebi_compounds(joinpath("data", "chebi", "chemical_data.tsv"), joinpath("data", "rhea", "rhea-compounds.txt"))
