module Lysine_Biosynthesis

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Lysine_Biosynthesis"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 34171, name = "4-hydroxy-tetrahydrodipicolinate synthase", isozymes = [[(4, "A0A1B1EEF1"),],], subsystem = subsystem,), #4.3.3.7
    (rhea_id = 17325, name = "2,3,4,5-tetrahydropyridine-2,6-dicarboxylate N-succinyltransferase", isozymes = [[(3, "A0A1B1EEE8"),],], subsystem = subsystem,), #2.3.1.117
    (rhea_id = 22608, name = "Succinyl-diaminopimelate desuccinylase", isozymes = [[(2, "A0A1B1EE51"),],], subsystem = subsystem,), #3.5.1.18
    (rhea_id = 15393, name = "Diaminopimelate epimerase", isozymes = [[(2, "A0A1B1EFY4"),],], subsystem = subsystem,), #5.1.1.7
    (rhea_id = 23676, name = "UDP-N-acetylmuramoyl-L-alanyl-D-glutamate--2,6-diaminopimelate ligase", isozymes = [[(1, "A0A1B1E9I1"),],], subsystem = subsystem,), #6.3.2.13
    (rhea_id = 28374, name = "UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine ligase", isozymes = [[(1, "A0A1B1E9K0"),],], subsystem = subsystem,), #6.3.2.10
    (rhea_id = 15101, name = "Diaminopimelate decarboxylase", isozymes = [[(4, "A0A1B1EGM8"),],], subsystem = subsystem,), #4.1.1.20
    #(rhea_id = 12929, name = "Homocitrate synthase", isozymes = [[(1, "A0A1B1EL49"),],], subsystem = subsystem,), #2.3.3.14
    #no information for 1.17.1.8
)

#=
These ECs were present but not included in rxns:
(rhea_id = 23776?, name = "Aspartokinase?", isozymes = [[(1, "A0A1B1EEX2"),], [(1, "A0A1B1EF8"),] ] ??, subsystem = subsystem,), #EC 2.7.2.4 not sure about all of it, also in glycine_thre.. 
(rhea_id = 24284, name = "Aspartate-semialdehyde dehydrogenase", isozymes = [[(2, "A0A1B1EDS3"),], [(2, "A0A1B1EDU6"),],], subsystem = subsystem,), #EC 1.2.1.11, also in glycine_thre.. 
(rhea_id = ?, name = "?", isozymes = [[(4, "A0A1B1E9N6"),],], subsystem = subsystem,), #EC 1.1.1.3 not sure about all of it, also in glycine_thre.. 

=#

end # module
