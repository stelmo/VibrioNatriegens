module Pyrimidine_Metabolism

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Pyrimidine Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 20013, name = "Aspartate transcarbamylase", isozymes = [[(1, "A0A1B1EF32"),],], subsystem = subsystem,), #2.1.3.2
    (rhea_id = 24296, name = "Dihydroorotase", isozymes = [[(2, "A0A1B1EJI4"),],], subsystem = subsystem,), #3.5.2.3
    (rhea_id = 30187, name = "Dihydroorotate dehydrogenase", isozymes = [[(1, "A0A1B1ECC1"),],], subsystem = subsystem,), #1.3.5.2
    (rhea_id = 10380, name = "Orotate phosphoribosyltransferase", isozymes = [[(2, "A0A1B1E8U0"),],], subsystem = subsystem,), #2.4.2.10
    (rhea_id = 11596, name = "OMP decarboxylase", isozymes = [[(2, "A0A1B1EDC2"),],], subsystem = subsystem,), #4.1.1.23
    (rhea_id = 24400, name = "Uridylate kinase", isozymes = [[(6, "A0A1B1EEC6"),], [(1, "B7UU64"),], [(1, "C0MP25"),], [(1, "C0MP26"),], [(1, "G9HWN4"),],], subsystem = subsystem,), #2.7.4.22
    (rhea_id = 26426, name = "Cytidine triphosphate synthetase", isozymes = [[(4, "A0A1B1EF87"),],], subsystem = subsystem,), #6.3.4.2

)

#=
These ECs were present but not included in rxns:

EC 6.3.5.5 - unclear subunit information
    (rhea_id = 18633, name = "Carbamoyl-phosphate synthase", isozymes = [[(),],], subsystem = subsystem,),

EC 2.7.4.6 - 2 rhea ids
    (rhea_id = , name = "Nucleoside diphosphate kinase", isozymes = [[(4, "A0A1B1E9Y0"),],], subsystem = subsystem,),

EC 3.6.1.9 - 2 rhea ids
    (rhea_id = , name = "Nucleoside triphosphate pyrophosphatase", isozymes = [[(1, "A0A1B1EF84"),],], subsystem = subsystem,),
=#

end # module
