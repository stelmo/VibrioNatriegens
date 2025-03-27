[ci-img]: https://github.com/stelmo/VibrioNatriegens.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/stelmo/VibrioNatriegens.jl/actions/workflows/ci.yml

# Genome-scale metabolic model of _Vibrio natriegens_

[![CI status][ci-img]][ci-url]

This package builds a fully manually reconstructed genome-scale metabolic model of the halophilic bacterium _Vibrio natriegens_. 
The model is composed of 1283 reactions, 950 metabolites, and 939 genes. 
The model focusses on the primary metabolism the organism, and includes enzyme and ribosomal constraints. 

## Model characterization
At a glance, the model consists of:

| Attribute | Number |
|-----------|-------|
| Metabolic reactions | 803 |
| Transport reactions | 327 |
| Exchange reactions | 153 |
| Metabolic reactions with GRRs | 752 |
| Transport reactions with GRRs | 88 |
| Blocked reactions | 0 |
| Deadend metabolites | 0 |

The model has the following reaction cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Metabolic reactions* | *803* |
| kegg.reaction | 736 (92%) |
| metacyc.reaction |  |
| ecocyc.reaction | ? |
| reactome.reaction | ? |
| rhea.reaction | 736 (92%) |
| bigg.reaction | 539 (67%) |
| ec.code |  |
| SBO | 754 (94 %) |

The model has the following metabolite cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total metabolites* | *950* |
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
| *Total genes* | 939 |
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
