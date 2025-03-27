# build the readme
using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
using SBMLFBCModels
import AbstractFBCModels as A
using Distributed

model = VibrioNatriegens.build_model()
VibrioNatriegens.print_reactions(model)
VibrioNatriegens.print_metabolites(model)

# write the model to the main directory in standardized formats
m = convert(JSONFBCModels.JSONFBCModel, model)
AbstractFBCModels.save(m, "vibrio_natriegens.json")

m = convert(SBMLFBCModels.SBMLFBCModel, model)
AbstractFBCModels.save(m, "vibrio_natriegens.sbml")

# if nprocs() == 1
#     addprocs(6, exeflags = `--project=$(Base.active_project())`)
# end

# @everywhere begin
#     using Gurobi, COBREXA, ConstraintTrees, AbstractFBCModels
#     import ConstraintTrees as C
#     import AbstractFBCModels as A
#     using VibrioNatriegens
# end

rids = A.reactions(model)

transport_rxns = [
    filter(startswith("PERM_"), rids)
    filter(startswith("SYM_"), rids)
    filter(startswith("ANTI_"), rids)
    filter(startswith("ABC_"), rids)
    filter(startswith("PTS_"), rids)
    filter(startswith("DF_"), rids)
]
exchange_rxns = filter(startswith("EX_"), rids)
metabolic_rxns = setdiff(rids, [transport_rxns; exchange_rxns])

metabolic_rxn_grrs = [
    rid for rid in metabolic_rxns if !isnothing(A.reaction_gene_association_dnf(model, rid))
]
transport_rxn_grrs = [
    rid for rid in transport_rxns if !isnothing(A.reaction_gene_association_dnf(model, rid))
]

blocked_rxns = []

deadend_metabolites = []

sbo_reactions =
    [rid for rid in metabolic_rxns if haskey(model.reactions[rid].annotations, "SBO")]
bigg_reactions = [
    rid for
    rid in metabolic_rxns if haskey(model.reactions[rid].annotations, "bigg.reaction")
]
kegg_reactions = [
    rid for
    rid in metabolic_rxns if haskey(model.reactions[rid].annotations, "kegg.reaction")
]
rhea_reactions = [
    rid for
    rid in metabolic_rxns if haskey(model.reactions[rid].annotations, "rhea.reaction")
]
ec_code =
    [rid for rid in metabolic_rxns if haskey(model.reactions[rid].annotations, "ec-code")]

mids = A.metabolites(model)

gids = A.genes(model)

readme = """
[ci-img]: https://github.com/stelmo/VibrioNatriegens.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/stelmo/VibrioNatriegens.jl/actions/workflows/ci.yml

# Genome-scale metabolic model of _Vibrio natriegens_

[![CI status][ci-img]][ci-url]

This package builds a fully manually reconstructed genome-scale metabolic model of the halophilic bacterium _Vibrio natriegens_. 
The model is composed of $(A.n_reactions(model)) reactions, $(A.n_metabolites(model)) metabolites, and $(A.n_genes(model)) genes. 
The model focusses on the primary metabolism the organism, and includes enzyme and ribosomal constraints. 

## Model characterization
At a glance, the model consists of:

| Attribute | Number |
|-----------|-------|
| Metabolic reactions | $(length(metabolic_rxns)) |
| Transport reactions | $(length(transport_rxns)) |
| Exchange reactions | $(length(exchange_rxns)) |
| Metabolic reactions with GRRs | $(length(metabolic_rxn_grrs)) |
| Transport reactions with GRRs | $(length(transport_rxn_grrs)) |
| Blocked reactions | $(length(blocked_rxns)) |
| Deadend metabolites | $(length(deadend_metabolites)) |

The model has the following reaction cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Metabolic reactions* | *$(length(metabolic_rxns))* |
| kegg.reaction | $(length(kegg_reactions)) ($( round(Int, length(kegg_reactions) / length(metabolic_rxns) * 100))%) |
| metacyc.reaction |  |
| ecocyc.reaction | ? |
| reactome.reaction | ? |
| rhea.reaction | $(length(rhea_reactions)) ($( round(Int, length(rhea_reactions) / length(metabolic_rxns) * 100))%) |
| bigg.reaction | $(length(bigg_reactions)) ($( round(Int, length(bigg_reactions) / length(metabolic_rxns) * 100))%) |
| ec.code |  |
| SBO | $(length(sbo_reactions)) ($( round(Int, length(sbo_reactions) / length(metabolic_rxns) * 100)) %) |

The model has the following metabolite cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total metabolites* | *$(length(mids))* |
| chebi.metabolite | ? |
| kegg.compound | ? |
| inchi | ? |
| inchikey | ? |
| smiles | ? |
| formula | ? |
| molarmass | ? |
| SBO | ? |


The model has the following gene cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total genes* | $(length(gids)) |
|  | ? |
| SBO | ? |


## Build instructions
For convenience, an SBML and JSON formatted model is available in the main directory of the repo.
If you would like to build the model from scratch (recommended, since some useful data gets lost converting to the commonly used model formats), then you must first install `VibrioNatriegens` like any Julia package.
```julia
] add VibrioNatriegens
```
Basic tests (MEMOTE-like) can be run by testing the package just like any other Julia package. 
These tests also run with CI every time the model is changed.
```julia

] test VibrioNatriegens
```
Thereafter, you only need to call a single function to build the entire model.
```julia
using VibrioNatriegens

model = VibrioNatriegens.build_model()
```
This model works well with the COBREXA package, and it can be used to simulate FBA, ec-FBA, and sRBA models.

#### Acknowledgements

`DifferentiableMetabolism.jl` was developed at Institute for Quantitative and
Theoretical Biology at Heinrich Heine University Düsseldorf
([qtb.hhu.de](https://www.qtb.hhu.de/en/)).

<img src="scripts/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="scripts/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto"> 

Note: this readme is automatically built using the script `scripts/readme.jl`.
"""

open("readme.md", "w") do io
    write(io, readme)
end
