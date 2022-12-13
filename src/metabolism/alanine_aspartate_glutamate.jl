module Ala_Asp_Glu

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Alanine Aspartate and Glutamate metabolism"

# (rhea_id = , name = , isozymes = [[(),],], subsystem = subsystem,),

rxns = (
    (rhea_id = 20013, name = "Aspartate carbamoyltransferase", isozymes = [[(1, "A0A1B1EF32"),],], subsystem = subsystem,), #EC 2.1.3.2
    (rhea_id = 15753, name = "Adenylosuccinate synthetase", isozymes = [[(2, "A0A1B1EFK9"),],[(2, "A0A1B1EHC6"),],], subsystem = "$subsystem, Purine Metabolism"), #6.3.4.4
    #(rhea_id = 10932, name = "Argininosuccinate synthase", isozymes = [[(3, "A0A1B1EFQ0"),],], subsystem = subsystem,), #6.3.4.5
    (rhea_id = 18405, name = "Alanine dehydrogenase", isozymes = [[(1, "A0A1B1EBB0"),],], subsystem = subsystem,), #1.4.1.1
    (rhea_id = 16601, name = "Aspartate ammonia-lyase", isozymes = [[(1, "A0A1B1EFP1"),],], subsystem = subsystem,), #4.3.1.1
    #(rhea_id = 24020, name = "Argininosuccinate lyase", isozymes = [[(1, "A0A1B1EG54"),],], subsystem = subsystem,), #4.3.2.1
    (rhea_id = 25876, name = "L-aspartate oxidase", isozymes = [[(1, "A0A1B1EF02"),],], subsystem = subsystem,), #1.4.3.16
    #(rhea_id = 15133, name = "Glutamate dehydrogenase", isozymes = [[(1, "A0A1B1ECC8"),],], subsystem = subsystem,), #1.4.1.2
    (rhea_id = 30235, name = "Bifunctional protein PutA", isozymes = [[(1, "A0A1B1EHY2"),],], subsystem = "$subsystem, Alanine Aspartate and Glutamate metabolism" ), #1.2.1.88
    (rhea_id = 15889, name = "Glutaminase", isozymes = [[(4, "A0A1B1EF16"),],], subsystem = subsystem,), #3.5.1.2
    (rhea_id = 13237, name = "Glutamine--fructose-6-phosphate aminotransferase", isozymes = [[(2, "A0A1B1E9C4"),],], subsystem = subsystem,), #2.6.1.16
    (rhea_id = 14905, name = "Amidophosphoribosyltransferase", isozymes = [[(1, "A0A1B1EDU7"),],], subsystem = "$subsystem, Purine Metabolism",), #2.4.2.14
)

#=
These ECs were present but not included in rxns:
(rhea_id = unknown, name = "UDP-3-O-acyl-GlcNAc deacetylase", isozymes = [[(1, "A0A1B1EDT5"),],], subsystem = subsystem,), #3.5.1.1
(rhea_id = #23920, #16853, name = "Adenylosuccinate lyase", isozymes = [[(1, "A0A1B1EBC6"),],], subsystem = "$subsystem, Purine Metabolism"), #4.3.2.2
(rhea_id = #16169, name = "Formate-dependent phosphoribosylglycinamide formyltransferase", isozymes = [[(?, "A0A1B1E8Q4"),],], subsystem = subsystem,), #6.3.1.2 unknown information
(rhea_id = 18633, name = "Carbamoyl-phosphate synthase large chain", isozymes = [[(1, "A0A1B1E9K5"),[(1, "A0A1B1E9L6"),],], subsystem = subsystem,), #6.3.5.5 strange subunit information
=#

end # module
