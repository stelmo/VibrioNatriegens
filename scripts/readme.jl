# build the readme
using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
using SBMLFBCModels
import AbstractFBCModels as A

model = VibrioNatriegens.build_model()

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

metabolic_rxn_grrs = [rid for rid in metabolic_rxns if !isnothing(A.reaction_gene_association_dnf(model, rid))] 
transport_rxn_grrs = [rid for rid in transport_rxns if !isnothing(A.reaction_gene_association_dnf(model, rid))]

blocked_rxns = []

mids = A.metabolites(model)

deadend_metabolites = []

gids = A.genes(model)

readme = """
# Vibrio natriegens

This package builds a genome-scale metabolic model of the halophilic bacterium Vibrio natriegens. 
The model is composed of $(A.n_reactions(model)) reactions, $(A.n_metabolites(model)) metabolites, and $(A.n_genes(model)) genes. 
The model focusses on the primary metabolism the organism, and includes enzyme and ribosomal constraints. 

## Model characterization
At a glance, the model consists of:

| Attribute | Value |
|-----------|-------|
| Metabolic reactions | $(length(metabolic_rxns)) |
| Transport reactions | $(length(transport_rxns)) |
| Exchange reactions | $(length(exchange_rxns)) |
| Metabolic reactions with GRRs | $(length(metabolic_rxn_grrs)) |
| Transport reactions with GRRs | $(length(transport_rxn_grrs)) |
| Blocked reactions | $(length(blocked_rxns)) |
| Deadend metabolites | $(length(deadend_metabolites)) |

The model has the following reaction cross-references (available under the `annotations` field):

| Attribute | Value |
|-----------|-------|
| Metabolic reactions | $(length(metabolic_rxns)) |
|-----------|-------|
| kegg.reaction | ? |
| metacyc.reaction | ? |
| ecocyc.reaction | ? |
| reactome.reaction | ? |
| rhea.reaction | ? |
| bigg.reaction | ? |
| ec.code | ? |

The model has the following metabolite cross-references (available under the `annotations` field):

| Attribute | Value |
|-----------|-------|
| Total metabolites | $(length(mids)) |
|-----------|-------|
| chebi.metabolite | ? |
| kegg.compound | ? |
| inchi | ? |
| inchikey | ? |
| smiles | ? |
| formula | ? |
| molarmass | ? |


The model has the following gene cross-references (available under the `annotations` field):

| Attribute | Value |
|-----------|-------|
| Total genes | $(length(gids)) |
|-----------|-------|
|  | ? |


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
"""

open("readme.md", "w") do io
    write(io, readme)
end
