module Arigine_Biosynthesis

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Arigine Biosynthesis"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 19597, name = "Arginine deiminase", isozymes = [[(1, "A0A1B1EF55"),],], subsystem = subsystem,), #3.5.3.6
    (rhea_id = 15941, name = "Acetylornithine deacetylase", isozymes = [[(2, "A0A1B1EFF3"),],], subsystem = subsystem,), #3.5.1.16
    (rhea_id = 19513, name = "Ornithine carbamoyltransferase", isozymes = [[(1, "A0A1B1EF60"),],], subsystem = subsystem,), #2.1.3.3
    (rhea_id = 14629, name = "Acetylglutamate kinase", isozymes = [[(1, "A0A1B1EFC5"),],], subsystem = subsystem,), #2.7.2.8
    (rhea_id = 21588, name = "N-acetyl-gamma-glutamyl-phosphate reductase", isozymes = [[(1, "A0A1B1EFC7"),],], subsystem = subsystem,), #1.2.1.38
    (rhea_id = 18049, name = "Acetylornithine aminotransferase", isozymes = [[(2, "A0A1B1EFJ8"),],], subsystem = subsystem,), #2.6.1.11
    (rhea_id = 24292, name = "Amino-acid acetyltransferase", isozymes = [[(1, "A0A1B1EEH5"),],], subsystem = subsystem,), #2.3.1.1
)

#=
These ECs were present but not included in rxns:

These ECs were already recorded in other subsystem files: EC6.3.4.5, EC4.3.2.1, EC3.5.1.2, EC6.3.1.2, EC1.4.1.2 and EC3.5.1.5
=#

end # module
