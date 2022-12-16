module Arg_Pro_Meta

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Arginine and Proline Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 10812 , name = "N-succinylglutamate 5-semialdehyde dehydrogenase", isozymes = [[(1, "A0A1B1EFI6"),],], subsystem = subsystem,), #1.2.1.71
    (rhea_id = 15169 , name = "Succinylglutamate desuccinylase", isozymes = [[(1, "A0A1B1EBU0"),],], subsystem = subsystem,), #3.5.1.96
    (rhea_id = 14877 , name = "Glutamate 5-kinase", isozymes = [[(1, "A0A1B1EA52"),],], subsystem = subsystem,), #2.7.2.11
    (rhea_id = 19541 , name = "Gamma-glutamyl phosphate reductase", isozymes = [[(1, "A0A1B1EA60"),],], subsystem = subsystem,), #1.2.1.41
    (rhea_id = 17641 , name = "Biosynthetic arginine decarboxylase", isozymes = [[(1, "A0A1B1ELE0"),],], subsystem = subsystem,), #4.1.1.19
    (rhea_id = 18037 , name = "Putative agmatine deiminase", isozymes = [[(1, "A0A1B1ECY6"),],], subsystem = subsystem,), #3.5.3.12
    (rhea_id = 23784 , name = "Proline dehydrogenase", isozymes = [[(1, "A0A1B1EHY2"),],], subsystem = subsystem,), #1.5.5.2 - part of the bifunctional protein PutA
)


#=
These ECs were present but not included in rxns:
    (rhea_id = 14105/14109 , name = "Pyrroline-5-carboxylate reductase", isozymes = [[(1, "A0A1B1EF33"),],], subsystem = subsystem,), #1.5.1.2 unknown rhea
    (rhea_id = 34099/34095 , name = "Carboxynorspermidine/carboxyspermidine decarboxylase", isozymes = [[(2, "A0A1B1ED48"),],], subsystem = subsystem,), #4.1.1.96 unknown rhea

These ECs were already recorded in other subsystem files: EC1.2.1.88
=#

end # module