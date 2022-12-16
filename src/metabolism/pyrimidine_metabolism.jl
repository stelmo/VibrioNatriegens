module Pyrimidine_Metabolism

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Pyrimidine Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 20013, name = "Aspartate transcarbamylase", isozymes = [[(1, "A0A1B1EF32"),],], subsystem = "$subsystem, Alanine, Aspartate and Glutamate Metabolism",), #2.1.3.2
    (rhea_id = 24296, name = "Dihydroorotase", isozymes = [[(2, "A0A1B1EJI4"),],], subsystem = subsystem,), #3.5.2.3
    (rhea_id = 30187, name = "Dihydroorotate dehydrogenase", isozymes = [[(1, "A0A1B1ECC1"),],], subsystem = subsystem,), #1.3.5.2
    (rhea_id = 10380, name = "Orotate phosphoribosyltransferase", isozymes = [[(2, "A0A1B1E8U0"),],], subsystem = subsystem,), #2.4.2.10
    (rhea_id = 11596, name = "OMP decarboxylase", isozymes = [[(2, "A0A1B1EDC2"),],], subsystem = subsystem,), #4.1.1.23
    (rhea_id = 24400, name = "Uridylate kinase", isozymes = [[(6, "A0A1B1EEC6"),], [(1, "B7UU64"),], [(1, "C0MP25"),], [(1, "C0MP26"),], [(1, "G9HWN4"),],], subsystem = subsystem,), #2.7.4.22
    (rhea_id = 26426, name = "Cytidine triphosphate synthetase", isozymes = [[(4, "A0A1B1EF87"),],], subsystem = subsystem,), #6.3.4.2
    (rhea_id = 23252, name = "Ribonucleoside-diphosphate reductase", isozymes = [[(1, "A0A1B1ED46"),], [(1, "A0A1B1ED65"),], [(1, "A0A5B9N9N9"),], [(1, "A0A5B9NCJ5"),],], subsystem = "$subsystem, Purine Metabolism",), #1.17.4.1
    (rhea_id = 13517, name = "Thymidylate kinase", isozymes = [[(1, "A0A1B1EDG3"),],], subsystem = subsystem,), #2.7.4.9
    (rhea_id = 12104, name = "Thymidylate synthase", isozymes = [[(2, "A0A1B1E9Q1"),],], subsystem = subsystem,), #2.1.1.45
    (rhea_id = 16825, name = "Uridine kinase", isozymes = [[(1, "A0A1B1EGH9"),],], subsystem = subsystem,), #2.7.1.48
    (rhea_id = 12484, name = "5'-nucleotidase SurE", isozymes = [[(1, "A0A1B1EEY0"),], [(1, "A0A1B1ELI9"),],], subsystem = "$subsystem, Purine Metabolism",), #3.1.3.5
    (rhea_id = 13017, name = "Uracil phosphoribosyltransferase", isozymes = [[(1, "A0A1B1EE91"),],], subsystem = subsystem,), #2.4.2.9
    (rhea_id = 36167, name = "5'-deoxynucleotidase", isozymes = [[(2, "A0A1B1EAF9"),],], subsystem = "$subsystem, Purine Metabolism",), #3.1.3.89
    (rhea_id = 19129, name = "Thymidine kinase", isozymes = [[(4, "A0A1B1EBF7"),],], subsystem = subsystem,), #2.7.1.21
    (rhea_id = 16069, name = "Cytidine deaminase", isozymes = [[(2, "A0A1B1EBR9"),],], subsystem = subsystem,), #3.5.4.5
    (rhea_id = 52540, name = "Pyrimidine nucleoside phosphorylase", isozymes = [[(1, "A0A1B1EHW8"),],], subsystem = subsystem,), #2.4.2.2
    (rhea_id = 24388, name = "Uridine phosphorylase", isozymes = [[(1, "A0A1B1EAY8"),],], subsystem = subsystem,), #2.4.2.3
    (rhea_id = 16037, name = "Thymidine phosphorylase", isozymes = [[(2, "A0A1B1EEL8"),],], subsystem = subsystem,), #2.4.2.4
)

#=
These ECs were present but not included in rxns:

EC 6.3.5.5 - unclear subunit information
    (rhea_id = 18633, name = "Carbamoyl-phosphate synthase", isozymes = [[(),],], subsystem = "$subsystem, Alanine, Aspartate and Glutamate Metabolism",),

EC 2.7.4.6 - 2 rhea ids
    (rhea_id = , name = "Nucleoside diphosphate kinase", isozymes = [[(4, "A0A1B1E9Y0"),],], subsystem = subsystem,),

EC 3.6.1.9 - 2 rhea ids
    (rhea_id = , name = "Nucleoside triphosphate pyrophosphatase", isozymes = [[(1, "A0A1B1EF84"),],], subsystem = "$subsystem, Purine Metabolism",),

EC 2.7.4.25 - 2 rhea ids
    (rhea_id = , name = "Cytidylate kinase", isozymes = [[(1, "A0A1B1EDN8"),],], subsystem = subsystem,),
=#

end # module
