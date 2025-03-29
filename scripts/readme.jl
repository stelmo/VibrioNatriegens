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

if nprocs() == 1
    addprocs(7, exeflags = `--project=$(Base.active_project())`)
end

@everywhere begin
    using Gurobi, COBREXA, ConstraintTrees, AbstractFBCModels
    import ConstraintTrees as C
    import AbstractFBCModels as A
    using VibrioNatriegens
end

include(joinpath("characterization_functions.jl"))

# reactions
rids = A.reactions(model)

transport_rxns = [
    filter(startswith("PERM_"), rids)
    filter(startswith("SYM_"), rids)
    filter(startswith("ANTI_"), rids)
    filter(startswith("ABC_"), rids)
    filter(startswith("PTS_"), rids)
]
exchange_rxns = filter(startswith("EX_"), rids)
diffusion_rxns = filter(startswith("DF_"), rids)
metabolic_rxns = setdiff(rids, [transport_rxns; exchange_rxns; diffusion_rxns])

metabolic_rxn_grrs = [
    rid for rid in metabolic_rxns if !isnothing(A.reaction_gene_association_dnf(model, rid))
]
transport_rxn_grrs = [
    rid for rid in transport_rxns if !isnothing(A.reaction_gene_association_dnf(model, rid))
]

blocked_rxns = blocked_reactions()

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

setdiff(A.reactions(model), rhea_reactions)

ec_reactions = [
    rid for
    rid in metabolic_rxns if haskey(model.reactions[rid].annotations, "eggnog.ec") ||
    haskey(model.reactions[rid].annotations, "rhea.ec") ||
    haskey(model.reactions[rid].annotations, "kegg.ec")
]

go_reactions =
    [rid for rid in metabolic_rxns if haskey(model.reactions[rid].annotations, "eggnog.go")]

metacyc_reactions = [
    rid for rid in metabolic_rxns if
    haskey(model.reactions[rid].annotations, "metacyc.reaction")
]

reactome_reactions = [
    rid for rid in metabolic_rxns if
    haskey(model.reactions[rid].annotations, "reactome.reaction")
]

seed_reactions = [
    rid for
    rid in metabolic_rxns if haskey(model.reactions[rid].annotations, "seed.reaction")
]

# metabolites
mids = A.metabolites(model)

formula_metabolites =
    [mid for mid in A.metabolites(model) if !isnothing(A.metabolite_formula(model, mid))]

charge_metabolites =
    [mid for mid in A.metabolites(model) if !isnothing(A.metabolite_charge(model, mid))]

sbo_metabolites = [
    mid for mid in A.metabolites(model) if haskey(model.metabolites[mid].annotations, "SBO")
]

inchi_metabolites = [
    mid for
    mid in A.metabolites(model) if haskey(model.metabolites[mid].annotations, "inchi")
]

inchikey_metabolites = [
    mid for mid in A.metabolites(model) if
    haskey(model.metabolites[mid].annotations, "inchikey")
]

smiles_metabolites = [
    mid for
    mid in A.metabolites(model) if haskey(model.metabolites[mid].annotations, "smiles")
]

chebi_metabolites = [
    mid for
    mid in A.metabolites(model) if haskey(model.metabolites[mid].annotations, "chebi")
]

kegg_metabolites = [
    mid for mid in A.metabolites(model) if
    haskey(model.metabolites[mid].annotations, "kegg.compound")
]

# genes
gids = A.genes(model)

sbo_genes = [gid for gid in A.genes(model) if haskey(model.genes[gid].annotations, "SBO")]
ncbi_genes = [
    gid for
    gid in A.genes(model) if haskey(model.genes[gid].annotations, "protein.accession")
]

readme = """
[ci-img]: https://github.com/stelmo/VibrioNatriegens/workflows/CI/badge.svg
[ci-url]: https://github.com/stelmo/VibrioNatriegens/actions/workflows/ci.yml

# A genome-scale metabolic model of _Vibrio natriegens_

This package builds a fully manually reconstructed genome-scale metabolic model of the halophilic bacterium:

<p align="center">
<i> 
Vibrio natriegens  
</i> [![CI status][ci-img]][ci-url]
</p> 

The model is composed of $(A.n_reactions(model)) reactions, $(A.n_metabolites(model)) metabolites, and $(A.n_genes(model)) genes. 
It focusses on the primary metabolism of _V. natriegens_, and includes data to facilitate the construction of enzyme and ribosomal constraints. 
A MEMOTE-like test suite is implemented in the `test` directory, and runs with CI. 
The primary name spaces for the model reactions and metabolites are [Rhea](https://www.rhea-db.org/) and [ChEBI](https://www.ebi.ac.uk/chebi/).
The RefSeq annotated genome was used to map genes to reactions can be downloaded from the NCBI accession [GCF_001456255.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001456255.1/
)

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

## Cross references

The model has the following reaction cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Metabolic reactions* | *$(length(metabolic_rxns))* |
| rhea.reaction | $(length(rhea_reactions)) ($( round(Int, length(rhea_reactions) / length(metabolic_rxns) * 100))%) |
| kegg.reaction | $(length(kegg_reactions)) ($( round(Int, length(kegg_reactions) / length(metabolic_rxns) * 100))%) |
| metacyc.reaction |  $(length(metacyc_reactions)) ($( round(Int, length(metacyc_reactions) / length(metabolic_rxns) * 100))%) |
| reactome.reaction | $(length(reactome_reactions)) ($( round(Int, length(reactome_reactions) / length(metabolic_rxns) * 100))%) |
| seed.reaction | $(length(seed_reactions)) ($( round(Int, length(seed_reactions) / length(metabolic_rxns) * 100))%) |
| eggnog.go | $(length(go_reactions)) ($( round(Int, length(go_reactions) / length(metabolic_rxns) * 100))%) |
| bigg.reaction | $(length(bigg_reactions)) ($( round(Int, length(bigg_reactions) / length(metabolic_rxns) * 100))%) |
| ec | $(length(ec_reactions)) ($( round(Int, length(kegg_reactions) / length(metabolic_rxns) * 100))%) |
| SBO | $(length(sbo_reactions)) ($( round(Int, length(sbo_reactions) / length(metabolic_rxns) * 100)) %) |

ec = `eggnog.ec` or `rhea.ec` or `kegg.ec`

The model has the following metabolite cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total metabolites* | *$(length(mids))* |
| chebi.metabolite | $(length(chebi_metabolites)) ($( round(Int, length(chebi_metabolites) / length(mids) * 100)) %) |
| kegg.compound | $(length(kegg_metabolites)) ($( round(Int, length(kegg_metabolites) / length(mids) * 100)) %) |
| inchi | $(length(inchi_metabolites)) ($( round(Int, length(inchi_metabolites) / length(mids) * 100)) %) |
| inchikey | $(length(inchikey_metabolites)) ($( round(Int, length(inchikey_metabolites) / length(mids) * 100)) %) |
| smiles | $(length(smiles_metabolites)) ($( round(Int, length(smiles_metabolites) / length(mids) * 100)) %) |
| formula | $(length(formula_metabolites)) ($( round(Int, length(formula_metabolites) / length(mids) * 100)) %) |
| charge | $(length(charge_metabolites)) ($( round(Int, length(charge_metabolites) / length(mids) * 100)) %) |
| SBO | $(length(sbo_metabolites)) ($( round(Int, length(sbo_metabolites) / length(mids) * 100)) %) |

The model has the following gene cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total genes* | $(length(gids)) |
| NCBI accession | $(length(ncbi_genes)) ($( round(Int, length(ncbi_genes) / length(gids) * 100)) %) |
| SBO | $(length(sbo_genes)) ($( round(Int, length(sbo_genes) / length(gids) * 100)) %) |

## Visualization
In the `maps` directory, [Escher compatible maps](https://escher.github.io/#/) of all the major subsystems of _V. natriegens_ may be found. 
Every metabolic reaction in the model is accounted for in these maps, facilitating a visual inspection of flux solutions etc.

## Build instructions

For convenience, an SBML and JSON formatted model is available in the main directory of the repo.

If you would like to build the model from scratch (recommended, since some useful data gets lost converting to the commonly used model formats), then you must first install `VibrioNatriegens` like any Julia package.
```julia
] add VibrioNatriegens
```
Note, the Rhea reactions used in the model are downloaded, and cached, from `https://www.rhea-db.org/` the first time the model is built, so it might take a few minutes at first.
Thereafter the building process should take less than 10 seconds.

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
The reactions, metabolites, and genes of `model` can be accessed like any field, e.g. `model.reactions`.

This model works well with the [COBREXA package](https://github.com/COBREXA/COBREXA.jl), and it can be used to simulate FBA, ec-FBA, and sRBA models.

## Known issues
Known issues in the model are marked as broken tests, mainly tested in the file `test/issues.jl` in the test directory.
Exceptions to this include substrate utilization tests (also listed as broken tests, but found in `test/basics.jl`):
- Formate
- Maltose
- Starch

These issues will be fixed in due course.

## Acknowledgements

`VibrioNatriegens` was developed at Institute for Quantitative and
Theoretical Biology at Heinrich Heine University Düsseldorf
([qtb.hhu.de](https://www.qtb.hhu.de/en/)).

<img src="scripts/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="scripts/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto"> 

Note: this readme is automatically built using the script `scripts/readme.jl`.
"""

open("readme.md", "w") do io
    write(io, readme)
end
