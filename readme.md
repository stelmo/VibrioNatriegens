[ci-img]: https://github.com/stelmo/VibrioNatriegens/workflows/CI/badge.svg
[ci-url]: https://github.com/stelmo/VibrioNatriegens/actions/workflows/ci.yml

# A genome-scale metabolic model of _Vibrio natriegens_

This package builds a fully manually reconstructed genome-scale metabolic model of the halophilic bacterium:

<p align="center">
<i> 
Vibrio natriegens
</i>
</p> 

The model is composed of 1388 reactions, 1020 metabolites, and 939 genes. 
It focusses on the primary metabolism the organism, and includes enzyme and ribosomal constraints. 

A MEMOTE-like test suite is implemented in the `test` directory, and runs with CI. 
Its current status is: [![CI status][ci-img]][ci-url]. 
The primary name spaces for the model reactions and metabolites are [Rhea](https://www.rhea-db.org/) and [ChEBI](https://www.ebi.ac.uk/chebi/). 

## Model characterization
At a glance, the model consists of:

| Attribute | Number |
|-----------|-------|
| Metabolic reactions | 803 |
| Transport reactions | 209 |
| Exchange reactions | 188 |
| Metabolic reactions with GRRs | 752 |
| Transport reactions with GRRs | 88 |
| Blocked reactions | 24 |

## Cross references

The model has the following reaction cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Metabolic reactions* | *803* |
| rhea.reaction | 801 (100%) |
| kegg.reaction | 776 (97%) |
| metacyc.reaction |  718 (89%) |
| reactome.reaction | 212 (26%) |
| seed.reaction | 770 (96%) |
| eggnog.go | 432 (54%) |
| bigg.reaction | 562 (70%) |
| ec | 786 (97%) |
| SBO | 803 (100 %) |

ec = `eggnog.ec` or `rhea.ec` or `kegg.ec`

The model has the following metabolite cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total metabolites* | *1020* |
| chebi.metabolite | 1020 (100 %) |
| kegg.compound | 348 (34 %) |
| inchi | 961 (94 %) |
| inchikey | 961 (94 %) |
| smiles | 961 (94 %) |
| formula | 1020 (100 %) |
| charge | 1020 (100 %) |
| SBO | 1020 (100 %) |

The model has the following gene cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total genes* | 939 |
| NCBI accession | 939 (100 %) |
| SBO | 939 (100 %) |

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
