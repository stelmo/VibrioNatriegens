# Genome-scale metabolic model of _Vibrio natriegens_

This package builds a fully manually reconstructed genome-scale metabolic model of the halophilic bacterium _Vibrio natriegens_. 
The model is composed of 1202 reactions, 887 metabolites, and 989 genes. 
The model focusses on the primary metabolism the organism, and includes enzyme and ribosomal constraints. 

## Model characterization
At a glance, the model consists of:

| Attribute | Number |
|-----------|-------|
| Metabolic reactions | 784 |
| Transport reactions | 297 |
| Exchange reactions | 121 |
| Metabolic reactions with GRRs | 751 |
| Transport reactions with GRRs | 125 |
| Blocked reactions | 0 |
| Deadend metabolites | 0 |

The model has the following reaction cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Metabolic reactions* | *784* |
| kegg.reaction | 737 (94%) |
| metacyc.reaction |  |
| ecocyc.reaction | ? |
| reactome.reaction | ? |
| rhea.reaction | 737 (94%) |
| bigg.reaction | 540 (69%) |
| ec.code |  |
| SBO | 753 (96 %) |

The model has the following metabolite cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total metabolites* | *887* |
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
| *Total genes* | 989 |
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

Note: this readme is automatically built using the script `scripts/readme.jl`.
