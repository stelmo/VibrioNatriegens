module Phe_Tyr_Try

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Phenylalanine, Tyrosine and Tryptophan Biosynthesis"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 14717, name = "Phospho-2-dehydro-3-deoxyheptonate aldolase", isozymes = [[(1, "A0A1B1E9S9"),], [(1, "A0A1B1ED47"),], [(1, "A0A1B1EGX0"),],], subsystem = subsystem,), #2.5.1.54
    (rhea_id = 21968, name = "3-dehydroquinate synthase", isozymes = [[(1, "A0A1B1EG44"),],], subsystem = subsystem,), #4.2.3.4
    (rhea_id = 21096, name = "3-dehydroquinate dehydratase", isozymes = [[(10, "A0A1B1EFQ5"),], [(10, "A0A1B1EGR9"),],], subsystem = subsystem,), #4.2.1.10
    (rhea_id = 17737, name = "Shikimate dehydrogenase (NADP(+))", isozymes = [[(2, "A0A1B1EGS5"),],], subsystem = subsystem,), #1.1.1.25
    (rhea_id = 13121, name = "Shikimate kinase", isozymes = [[(1, "A0A1B1EFP2"),],], subsystem = subsystem,), #2.7.1.71
    (rhea_id = 21256, name = "3-phosphoshikimate 1-carboxyvinyltransferase", isozymes = [[(1, "A0A1B1EBS4"),],], subsystem = subsystem,), #2.5.1.19
    (rhea_id = 21020, name = "Chorismate synthase", isozymes = [[(4, "A0A1B1EDV6"),],], subsystem = subsystem,), #4.2.3.5
    (rhea_id = 21648, name = "Prephenate dehydratase", isozymes = [[(1, "A0A1B1E9S4"),],], subsystem = subsystem,), #4.2.1.51 - this and 5.4.99.5 form a bifunctional complex
    (rhea_id = 13897, name = "Chorismate mutase", isozymes = [[(1, "A0A1B1E9S4"),],], subsystem = subsystem,), #5.4.99.5 - this and 4.2.1.51 form a bifunctional complex
    (rhea_id = 21732, name = "Anthranilate synthase component 1", isozymes = [[(1, "A0A1B1ED69"),],], subsystem = subsystem,), #4.1.3.27
    (rhea_id = 11768, name = "Anthranilate phosphoribosyltransferase", isozymes = [[(2, "A0A1B1EDG9"),],], subsystem = subsystem,), #2.4.2.18
    (rhea_id = 21540, name = "N-(5'-phosphoribosyl)anthranilate isomerase", isozymes = [[(1, "A0A1B1ED64"),],], subsystem = subsystem,), #5.3.1.24 - this and 4.1.1.48 are 1 multifunctional fusion protein
    (rhea_id = 23476, name = "Indole-3-glycerol phosphate synthase", isozymes = [[(1, "A0A1B1ED64"),],], subsystem = subsystem,), #4.1.1.48 - this and 5.3.1.24 are 1 multifunctional fusion protein
)

#=
These ECs were already recorded in other subsystem files: EC4.2.1.20, EC2.6.1.9 and EC1.14.16.1
=#

end # module
