module Glycolysis_Gluconeogenesis

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Glycolysis/Gluconeogenesis"

# (rhea_id = , name = , isozymes = [[(),],], subsystem = subsystem,),
 
rxns = (
    (rhea_id = 23536, name = "Phosphoglucomutase", isozymes = [[(1, "A0A1B1EB72"),], [(1, "A0A1B1EER5"),]], subsystem = subsystem,),
    (rhea_id = 11816, name = "Glucose-6-phosphate isomerase", isozymes = [[(1, "A0A1B1EFC4"),],], subsystem = subsystem,),
    (rhea_id = 11064, name = "Fructose-1,6-bisphosphatase", isozymes = [[(4, "A0A1B1ELE3"),], [(4, "A0A1B1ELE3"),],], subsystem = subsystem,),
)

end # module
