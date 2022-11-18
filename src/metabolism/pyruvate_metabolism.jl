module Pyruvate_Metabolism

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Pyruvate Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 21864, name = "Hydroxyacylglutathione hydrolase", isozymes = [[(1, "A0A1B1EEA6"),],], subsystem = subsystem,), #3.1.2.6
    (rhea_id = 19069, name = "Lactoylglutathione lyase", isozymes = [[(1, "A0A1B1EDM9"),],], subsystem = subsystem,), #4.4.1.5
    (rhea_id = 51468, name = "D-lactate dehydrogenase", isozymes = [[(1, "A0A1B1EKM0"),],], subsystem = subsystem,), #1.1.5.12
    (rhea_id = 18157, name = "Pyruvate kinase", isozymes = [[(1, "A0A1B1E9J8"),], [(1, "A0A1B1EDE7"),],], subsystem = subsystem,), #2.7.1.40
    (rhea_id = 11364, name = "PEP synthase", isozymes = [[(1, "A0A1B1EHZ4"),], [(1, "A0A1B1EK98"),],], subsystem = subsystem,), #2.7.9.2
    (rhea_id = 28370, name = "Phosphoenolpyruvate carboxylase", isozymes = [[(4, "A0A1B1EFF8"),],], subsystem = subsystem,), #4.1.1.31
    (rhea_id = 18617, name = "PEP carboxykinase", isozymes = [[(1, "A0A1B1E8K2"),],], subsystem = subsystem,), #4.1.1.49
    (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #4.1.1.12
)

#=
These ECs were present but not included in rxns:

EC 7.2.4.2 - includes a heterotrimer
    (rhea_id = 57724, name = "Oxaloacetate decarboxylase", isozymes = [[(),],], subsystem = subsystem,),
    A0A1B1EEW0  SUBUNIT: Heterotrimer of an alpha, a beta and a gamma subunit. {ECO:0000256|ARBA:ARBA00011869, ECO:0000256|HAMAP-Rule:MF_00404}.
    A0A1B1EEU7  SUBUNIT: Heterotrimer of an alpha, a beta and a gamma subunit. {ECO:0000256|ARBA:ARBA00011869, ECO:0000256|PIRNR:PIRNR015658}.

=#

end # module
