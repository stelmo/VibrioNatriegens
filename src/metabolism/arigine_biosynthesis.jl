module Arigine_Biosynthesis

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Arigine_biosynthesis"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 10932, name = "Argininosuccinate synthase", isozymes = [[(4, "A0A1B1EFQ0"),],], subsystem = subsystem,), #6.3.4.5
    (rhea_id = 24020, name = "Imidazole glycerol phosphate synthase subunit HisH", isozymes = [[(1, "A0A1B1EG54"),],], subsystem = subsystem,), #4.3.2.1
    (rhea_id = 19597, name = "Arginine deiminase", isozymes = [[(1, "A0A1B1EF55"),],], subsystem = subsystem,), #3.5.3.6
    (rhea_id = 15941, name = "Acetylornithine deacetylase", isozymes = [[(2, "A0A1B1EFF3"),],], subsystem = subsystem,), #3.5.1.16
    (rhea_id = 19513, name = "Ornithine carbamoyltransferase", isozymes = [[(1, "A0A1B1EF60"),],], subsystem = subsystem,), #2.1.3.3
    (rhea_id = 15133, name = "Glutamate dehydrogenase", isozymes = [[(1, "A0A1B1ECC8"),],], subsystem = subsystem,), #1.4.1.2
    (rhea_id = 14629, name = "Acetylglutamate kinase", isozymes = [[(1, "A0A1B1EFC5"),],], subsystem = subsystem,), #2.7.2.8
    (rhea_id = 21588, name = "N-acetyl-gamma-glutamyl-phosphate reductase", isozymes = [[(1, "A0A1B1EFC7"),],], subsystem = subsystem,), #1.2.1.38
    (rhea_id = 18049, name = "Acetylornithine aminotransferase", isozymes = [[(2, "A0A1B1EFJ8"),],], subsystem = subsystem,), #2.6.1.11
)

#=
These ECs were present but not included in rxns:
(rhea_id = 16169, name = "Formate-dependent phosphoribosylglycinamide formyltransferase", isozymes = [[(?, "A0A1B1E8Q4"),],], subsystem = subsystem,), #EC 6.3.1.2 unknown subunit information
(rhea_id = 15889, name = "Imidazole glycerol phosphate synthase subunit HisH", isozymes = [[(4, "A0A1B1EF16"),],], subsystem = subsystem,), #3.5.1.2 unknown information
no information to 2.3.1.1
(rhea_id = 20557, name = Urease subunit alpha isozymes = [[(?),],], subsystem = subsystem,), #EC 3.5.1.5 Heterotrimer


=#

end # module
