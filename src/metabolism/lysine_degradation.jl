module Lysine_Degradation

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Lysine Degradation"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    
    
)

#=
These ECs were present but not included in rxns:

EC 4.2.1.17 - missing rhea id and Heterotetramer
    (rhea_id = , name = "Fatty acid oxidation complex subunit alpha", isozymes = [[(),],], subsystem = "$subsystem, Tryptophan Metabolism, Phenylalanine Metabolism, Valine, Leucine and Isoleucine Biosynthesis",),

EC 1.1.1.35 unsure about rhea id and Heterotetramer
    (rhea_id = 21760?, name = "Fatty acid oxidation complex subunit alpha", isozymes = [[(?),],], subsystem = "$subsystem, Tryptophan Metabolism, Valine, Leucine and Isoleucine Biosynthesis",),

EC 3.4.-.-
    (rhea_id = , name = "Putative beta-barrel assembly-enhancing protease", isozymes = [[(1, "A0A1B1EEI3"),], [(1, "A0A1B1ELG0"),],], subsystem = subsystem,), 

EC 1.2.4.2 - missing rhea id, possibly 12188
    (rhea_id = , name = "Oxoglutarate dehydrogenase", isozymes = [[(),],], subsystem = "$subsystem, Citrate Cycle, Tryptophan Metabolism",),
=#

end # module
