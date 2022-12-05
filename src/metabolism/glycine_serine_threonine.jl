module Gly_Ser_Thr

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Glycine, Serine and Threonine Metabolism"

# (rhea_id = , name = , isozymes = [[(),],], subsystem = subsystem,),
 
rxns = (
    (rhea_id = 15901, name = "2,3-bisphosphoglycerate-independent phosphoglycerate mutase", isozymes = [[(1, "A0A1B1EFI9"),],], subsystem = subsystem,), #EC 5.4.2.12
    (rhea_id = 12641, name = "D-3-phosphoglycerate dehydrogenase", isozymes = [[(1, "A0A1B1EF12"),],], subsystem = subsystem,), #EC 1.1.1.95 cen
    (rhea_id = 19169, name = "L-serine dehydratase", isozymes = [[(1, "A0A1B1ED13"),], [(1, "A0A1B1EHL2"),], [(1, "A0A1B1EHL2"),],], subsystem = subsystem,), #EC 4.3.1.17
    (rhea_id = 13977, name = "Probable D-serine dehydratase", isozymes = [[(1, "A0A1B1EBV8"),],], subsystem = subsystem,), #EC 4.3.1.18
    (rhea_id = 22108, name = "L-threonine dehydratase", isozymes = [[(4, "A0A1B1EG77"),],], subsystem = subsystem,), #EC 4.3.1.19
    (rhea_id = 15481, name = "Serine hydroxymethyltransferase", isozymes = [[(2, "A0A1B1EAE3"),], [(2, "A0A1B1EJR1"),], [(2, "A0A1B1EJZ7"),],], subsystem = subsystem,), #EC 2.1.2.1
    (rhea_id = 22108, name = "Aminomethyltransferase", isozymes = [[(1, "A0A1B1EK07"),],], subsystem = subsystem,), #EC 2.1.2.10
    (rhea_id = 15045, name = "Dihydrolipoyl dehydrogenase", isozymes = [[(1, "A0A1B1EET9"),],], subsystem = subsystem,), #EC 1.8.1.4
    (rhea_id = 17433, name = "Oxygen-dependent choline dehydrogenase", isozymes = [[(1, "A0A1B1ELG7"),],], subsystem = subsystem,), #EC 1.1.99.1
    (rhea_id = 20736, name = "2-amino-3-ketobutyrate coenzyme A ligase", isozymes = [[(2, "A0A1B1EID3"),],], subsystem = subsystem,), #EC 2.3.1.29
    (rhea_id = 13161, name = "L-threonine 3-dehydrogenase", isozymes = [[(4, "A0A1B1EIG6"),],], subsystem = subsystem,), #EC 1.1.1.103
    (rhea_id = 10840, name = "Threonine synthase", isozymes = [[(1, "A0A1B1EAB0"),],], subsystem = subsystem,), #EC 4.2.3.1
    (rhea_id = 13985, name = "Homoserine kinase", isozymes = [[(1, "A0A1B1E9L7"),],], subsystem = subsystem,), #EC 2.7.1.39
    (rhea_id = 24284, name = "Aspartate-semialdehyde dehydrogenase", isozymes = [[(2, "A0A1B1EDS3"),], [(2, "A0A1B1EDU6"),],], subsystem = subsystem,), #EC 1.2.1.11
    (rhea_id = 11160, name = "Diaminobutyrate--2-oxoglutarate transaminase", isozymes = [[(1, "A0A1B1ECM4"),],], subsystem = subsystem,), #EC 2.6.1.76
    (rhea_id = 16901, name = "L-2,4-diaminobutyric acid acetyltransferase", isozymes = [[(1, "A0A1B1ECJ3"),],], subsystem = subsystem,), #EC 2.3.1.178
    (rhea_id = 17281, name = "L-ectoine synthase", isozymes = [[(1, "A0A1B1ECL0"),], [(1, "A0A1B1EIG0"),],], subsystem = subsystem,), #EC 4.2.1.108



)



#=
These ECs were present but not included in rxns:
(rhea_id = 16573/14329, name = "Phosphoserine aminotransferase", isozymes = [[(2, "A0A1B1ECA5"),],], subsystem = subsystem,), #EC 2.6.1.52
(rhea_id = 24873/21208, name = "Phosphoserine phosphatase", isozymes = [[(1, "A0A1B1EEW1"),],], subsystem = subsystem,), #EC 3.1.3.3
(rhea_id = 10532, name = "Tryptophan synthase beta chai", isozymes = [[(4, "A0A1B1ED60"),], [(4, "A0A1B1ED72"),], [(4, "A0A1B1EJ76"),],], subsystem = subsystem,), #EC 4.2.1.20 unknown subunit information 
(rhea_id = 24304, name = "Glycine dehydrogenase (decarboxylating)", isozymes = [[(?, "A0A1B1EK36"),],], subsystem = subsystem,), #EC 1.4.4.2 unkown subunit information
(rhea_id = 15305, name = "Betaine aldehyde dehydrogenase", isozymes = [[(?, "A0A1B1EL55"),],], subsystem = subsystem,), #EC 1.2.1.8 unkown subunit information
(rhea_id = ?, name = "?", isozymes = [[(4, "A0A1B1E9N6"),],], subsystem = subsystem,), #EC 1.1.1.3 not sure about all of it 
(rhea_id = 23776?, name = "Aspartokinase?", isozymes = [[(1, "A0A1B1EEX2"),], [(1, "A0A1B1EF8"),] ] ??, subsystem = subsystem,), #EC 2.7.2.4 not sure about all of it 

=#

end # module
