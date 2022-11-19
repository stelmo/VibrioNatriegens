module Citrate_Cycle

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Citrate Cycle"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 17017, name = "pyruvate dehydrogenase E2", isozymes = [[(1, "A0A1B1EES8"),],], subsystem = subsystem,), #2.3.1.12
    (rhea_id = 15045, name = "Dihydrolipoyl dehydrogenase", isozymes = [[(1, "A0A1B1EET9"),],], subsystem = subsystem,), #1.8.1.4
    (rhea_id = 18617, name = "Phosphoenolpyruvate carboxykinase", isozymes = [[(1, "A0A1B1E8K2"),],], subsystem = subsystem,), #4.1.1.49
    (rhea_id = 19629, name = "Isocitrate dehydrogenase", isozymes = [[(1, "A0A1B1EB33"),],], subsystem = subsystem,), #1.1.1.42
    (rhea_id = 15213, name = "2-oxoglutarate dehydrogenase E2", isozymes = [[(1, "A0A1B1EAS6"),],], subsystem = subsystem,), #2.3.1.61
    (rhea_id = 40523, name = "Succinate dehydrogenase", isozymes = [[(1, "A0A1B1EAJ6"),], [(1, "A0A1B1EAM6"),],], subsystem = subsystem,), #1.3.5.1
    (rhea_id = 12460, name = "Fumarate hydratase", isozymes = [[(4, "A0A1B1EFM8"),], [(2, "A0A1B1ED17"),],], subsystem = subsystem,), #4.2.1.2
    (rhea_id = 46012, name = "Malate dehydrogenase [quinone]", isozymes = [[(1, "A0A1B1EBX5"),],], subsystem = subsystem,), #1.1.5.4
    (rhea_id = 21432, name = "Malate dehydrogenase", isozymes = [[(2, "A0A1B1E9A0"),],], subsystem = subsystem,), #1.1.1.37
)

#=
These ECs were present but not included in rxns:

EC 1.2.4.1 - contains a heterodimer
    (rhea_id = 19189, name = "Pyruvate dehydrogenase E1", isozymes = [[(2, "A0A1B1EFH4"),], [(1, "A0A1B1EK05"),], [(1, "A0A1B1EK16"),],], subsystem = subsystem,),

EC 4.2.1.3 - 2 versions, one is associated with 2 EC numbers
    (rhea_id = 10336, name = "Aconitate hydratase", isozymes = [[(1, "A0A1B1EES2"),], [(1, "A0A1B1ECG9"),],], subsystem = subsystem,),

EC 1.2.4.2 - missing rhea_i_d, possibly 12188
    (rhea_id = , name = "Oxoglutarate dehydrogenase", isozymes = [[(1, "A0A1B1EB81"),],], subsystem = subsystem,),

EC 6.2.1.5 - contains a heterotetramer
    (rhea_id = 17661, name = "Succinyl-CoA synthetase subunit beta", isozymes = [[(4, "A0A1B1EAJ2"),], [(4, "A0A1B1EAK1"),],], subsystem = subsystem,),
=#

end # module
