module Histidine_Meta

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Histidine_Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 18473, name = "ATP phosphoribosyltransferase", isozymes = [[(1, "A0A1B1EBD2"),],], subsystem = subsystem,), #2.4.2.17
    (rhea_id = 22828, name = "Histidine biosynthesis bifunctional protein HisIE", isozymes = [[(1, "A0A1B1EBE1"),],], subsystem = subsystem,), #3.6.1.31?
    (rhea_id = 20049, name = "Histidine biosynthesis bifunctional protein HisIE", isozymes = [[(1, "A0A1B1EBE1"),],], subsystem = subsystem,), #3.5.4.19?
    (rhea_id = 15469, name = "1-(5-phosphoribosyl)-5-[(5-phosphoribosylamino)methylideneamino] imidazole-4-carboxamide isomerase", isozymes = [[(1, "A0A1B1EBE2"),],], subsystem = subsystem,), #5.3.1.16
    (rhea_id = 14465, name = "Histidine biosynthesis bifunctional protein HisB", isozymes = [[(1, "A0A1B1EBE7"),],], subsystem = subsystem,), #4.2.1.19 ?
    (rhea_id = 11040, name = "Histidine biosynthesis bifunctional protein HisB", isozymes = [[(1, "A0A1B1EBE7"),],], subsystem = subsystem,), #3.1.3.15 ?
    (rhea_id = 20641, name = "Histidinol dehydrogenase", isozymes = [[(1, "A0A1B1EBD6"),],], subsystem = subsystem,), #1.1.1.23
    (rhea_id = 21232, name = "Histidine ammonia-lyase", isozymes = [[(1, "A0A1B1EHA2"),],], subsystem = subsystem,), #4.3.1.3
    (rhea_id = 13101, name = "Urocanate hydratase", isozymes = [[(1, "A0A1B1EH77"),],], subsystem = subsystem,), #4.2.1.49
    (rhea_id = 23660, name = "Imidazolonepropionase", isozymes = [[(1, "A0A1B1EBS0"),],], subsystem = subsystem,), #3.5.2.7
    (rhea_id = 22492, name = "Formimidoylglutamase", isozymes = [[(1, "A0A1B1EBP9"),],], subsystem = subsystem,), #3.5.3.8

)


#=
These ECs were present but not included in rxns:
#strange information for 4.3.2.10
strange information for 2.1.1.-
=#

end # module