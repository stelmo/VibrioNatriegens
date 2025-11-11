[ci-img]: https://github.com/stelmo/VibrioNatriegens/workflows/CI/badge.svg
[ci-url]: https://github.com/stelmo/VibrioNatriegens/actions/workflows/ci.yml

# A genome-scale metabolic model of _Vibrio natriegens_

This package builds a fully manually reconstructed genome-scale metabolic model of the halophilic bacterium:

<p align="center">
<i> 
Vibrio natriegens  
</i> 
</p> 

The model is composed of 1421 reactions, 1022 metabolites, and 995 genes. 
It focusses on the primary metabolism of _V. natriegens_, and includes data to facilitate the construction of enzyme and ribosomal constraints. 
A MEMOTE-like test suite is implemented in the `test` directory, and runs with CI. 
The primary name spaces for the model reactions and metabolites are [Rhea](https://www.rhea-db.org/) and [ChEBI](https://www.ebi.ac.uk/chebi/).
The RefSeq annotated genome maps genes to reactions, and can be downloaded with the NCBI accession [GCF_001456255.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001456255.1/
)

## Model characterization
At a glance, the model consists of:

| Attribute | Value |
|-----------|-------|
| Metabolic reactions | 799 |
| Transport reactions | 242 |
| Exchange reactions | 189 |
| Metabolic reactions with GRRs | 750 |
| Transport reactions with GRRs | 122 |
| Blocked reactions | 36 |
| Test suite | [![CI status][ci-img]][ci-url] |

## Cross references

The model has the following reaction cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Metabolic reactions* | *799* |
| rhea.reaction | 799 (100%) |
| kegg.reaction | 774 (97%) |
| metacyc.reaction |  716 (90%) |
| reactome.reaction | 211 (26%) |
| seed.reaction | 661 (83%) |
| bigg.reaction | 601 (75%) |
| ec | 771 (97%) |
| SBO | 799 (100 %) |

Note: ec = `eggnog.ec` or `rhea.ec` or `kegg.ec`.


The model has the following metabolite cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total metabolites* | *1022* |
| chebi.metabolite | 1022 (100 %) |
| kegg.compound | 349 (34 %) |
| inchi | 963 (94 %) |
| inchikey | 963 (94 %) |
| smiles | 963 (94 %) |
| formula | 1022 (100 %) |
| charge | 1022 (100 %) |
| SBO | 1022 (100 %) |

The model has the following gene cross-references (available under the `annotations` field):

| Attribute | Number (fraction %) |
|-----------|-------|
| *Total genes* | 995 |
| NCBI accession | 995 (100 %) |
| SBO | 995 (100 %) |

## Visualization
In the `maps` directory, [Escher compatible maps](https://escher.github.io/#/) of all the major subsystems of _V. natriegens_ may be found. 
Every metabolic reaction in the model is accounted for in these maps, facilitating a visual inspection of flux solutions etc.

## Build instructions

For convenience, an SBML and JSON formatted model is available in the main directory of the repo.

If you would like to build the model from scratch (recommended, since some useful data gets lost converting to the commonly used model formats), then you must first install `VibrioNatriegens` like any Julia package.
```julia
] add https://github.com/stelmo/VibrioNatriegens.git
```
Note, the Rhea reactions used in the model are downloaded, and cached, from `https://www.rhea-db.org/` the first time the model is built, so it might take a few minutes at first.
Thereafter the building process should take less than 10 seconds.

Basic tests (MEMOTE-like) can be run by testing the package just like any other Julia package. 
These including mass, charge and stoichiometric consistency checks. 
Additionally, metabolite utilization, secretion, growth rates, and ATP yields are checked against measured values.
These tests run with CI every time the model is changed, but can be manually triggered with:
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

Currently the model does not simulate anaerobic growth accurately, as the ATP yield is only 3 ATP/glucose. 
This results in a low growth rate cf. experimental measurements. 
It is likely that an unknown anaerobic energy generating mechanism exists, but this is not captured in the model currently.

Additionally, the GAM and NGAM are set using experimental measurements, but they are relatively high cf. other metabolic reconstructions.
It is unclear if these are correct estimates, or not. Unfortunately most other models do not measure these values, but rather use default values supplied by some reconstruction tool.

These issues will be fixed in due course.

## Acknowledgements

`VibrioNatriegens` was developed at Institute for Quantitative and
Theoretical Biology at Heinrich Heine University Düsseldorf
([qtb.hhu.de](https://www.qtb.hhu.de/en/)).

<img src="scripts/assets/hhu.svg" alt="HHU logo" height="64px" style="height:64px; width:auto">   <img src="scripts/assets/qtb.svg" alt="QTB logo" height="64px" style="height:64px; width:auto"> 

Note: this readme is automatically built using the script `scripts/readme.jl`.
