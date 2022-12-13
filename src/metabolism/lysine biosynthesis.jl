module Lysine_Biosynthesis

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Lysine Biosynthesis"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 34171, name = "4-hydroxy-tetrahydrodipicolinate synthase", isozymes = [[(4, "A0A1B1EEF1"),],], subsystem = subsystem,), #4.3.3.7
    (rhea_id = 17325, name = "2,3,4,5-tetrahydropyridine-2,6-dicarboxylate N-succinyltransferase", isozymes = [[(3, "A0A1B1EEE8"),],], subsystem = subsystem,), #2.3.1.117
    (rhea_id = 22608, name = "Succinyl-diaminopimelate desuccinylase", isozymes = [[(2, "A0A1B1EE51"),],], subsystem = subsystem,), #3.5.1.18
    (rhea_id = 15393, name = "Diaminopimelate epimerase", isozymes = [[(2, "A0A1B1EFY4"),],], subsystem = subsystem,), #5.1.1.7
    (rhea_id = 23676, name = "UDP-N-acetylmuramoyl-L-alanyl-D-glutamate--2,6-diaminopimelate ligase", isozymes = [[(1, "A0A1B1E9I1"),],], subsystem = subsystem,), #6.3.2.13
    (rhea_id = 28374, name = "UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine ligase", isozymes = [[(1, "A0A1B1E9K0"),],], subsystem = subsystem,), #6.3.2.10
    (rhea_id = 15101, name = "Diaminopimelate decarboxylase", isozymes = [[(4, "A0A1B1EGM8"),],], subsystem = subsystem,), #4.1.1.20 
)

#=
These ECs were present but not included in rxns:
EC 1.17.1.8 - 2 rhea ids: RHEA:35323 and RHEA:35331
    (rhea_id = , name = "4-hydroxy-tetrahydrodipicolinate reductase", isozymes = [[(4, "A0A1B1E9L4"),],], subsystem = subsystem,),

These ECs were already recorded in other subsystem files: EC1.2.1.11, EC2.3.3.14, EC1.1.1.3, EC2.7.2.4
=#

end # module
