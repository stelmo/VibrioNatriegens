module Glycolysis_Gluconeogenesis

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Glycolysis/Gluconeogenesis"

# (rhea_id = , name = , isozymes = [[(),],], subsystem = subsystem,),
 
rxns = (
    (rhea_id = 23536, name = "Phosphoglucomutase", isozymes = [[(1, "A0A1B1EB72"),], [(1, "A0A1B1EER5"),]], subsystem = subsystem,), #EC 5.4.2.2
    (rhea_id = 11816, name = "Glucose-6-phosphate isomerase", isozymes = [ [(1, "A0A1B1EFC4"),], ], subsystem = subsystem,), #EC 5.3.1.9
    (rhea_id = 11064, name = "Fructose-1,6-bisphosphatase", isozymes = [[(4, "A0A1B1ELE3"),], [(4, "A0A1B1E9A1"),],], subsystem = subsystem,), #EC 3.1.3.11
    (rhea_id = 16109, name = "ATP-dependent 6-phosphofructokinase", isozymes = [[(4, "A0A1B1EFN6"),],], subsystem = subsystem,), #EC 2.7.1.11
    (rhea_id = 18585, name = "Triosephosphate isomerase", isozymes = [[(2, "A0A1B1E928"),],], subsystem = subsystem,), #EC 5.3.1.1
    (rhea_id = 14801, name = "Phosphoglycerate kinase", isozymes = [[(1, "A0A1B1EF14"),],], subsystem = subsystem,), #EC 2.7.2.3
)

#=
These ECs were present but not included in rxns:

2.7.1.199 - PTS transporter

=#

end # module
