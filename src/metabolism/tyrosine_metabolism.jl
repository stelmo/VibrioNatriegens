module Tyr_Metabolism

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Tyrosine Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 23744, name = "Histidinol-phosphate transaminase", isozymes = [[(2, "A0A1B1EBD9"),], [(2, "A0A1B1ELP6"),],], subsystem = "$subsystem, Histidine Metabolism, Phenylalanine Metabolism, Phenylalanine, Tyrosine and Tryptophan Biosynthesis",), #2.6.1.9
)

#=
These ECs were present but not included in rxns:

EC 2.3.1.- missing rhea id, name is possibly also Dihydrolipoamide acetyltransferase
    (rhea_id = , name = "Acetyltransferase", isozymes = [[(1, "A0A1B1EJX5"),], [(1, "A0A1B1EC85"),], [(1, "A0A1B1EJZ6"),],], subsystem = "$subsystem, Histidine Metabolism",),

EC 2.1.1.- confusing
    (rhea_id = , name = "", isozymes = [[(),],], subsystem = "$subsystem, Histidine Metabolism, Tryptophan Metabolism",),
=#

end # module
