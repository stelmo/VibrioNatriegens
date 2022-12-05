module Phe_Metabolism

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Phenylalanine Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 20273, name = "Phenylalanine-4-hydroxylase", isozymes = [[(1, "A0A1B1EJH6"),],], subsystem = "$subsystem, Phenylalanine, Tyrosine and Tryptophan Biosynthesis",), #1.14.16.1
    (rhea_id = 19481, name = "3-oxoadipyl-CoA thiolase", isozymes = [[(1, "A0A1B1EGV9"),],], subsystem = subsystem,), #2.3.1.174
    (rhea_id = 22624, name = "4-hydroxy-2-oxovalerate aldolase", isozymes = [[(1, "A0A1B1EC96"),],], subsystem = subsystem,), #4.1.3.39
)

#=
These ECs were present but not included in rxns:

EC 1.11.1.21 - confusing subunit info and 2 rhea ids: RHEA:20309 and RHEA:30275
    (rhea_id = , name = "Catalase-peroxidase", isozymes = [[(, "A0A1B1EJW7"),], [(, "A0A1B1EI40"),], [(, "A0A1B1EKV8"),],], subsystem = "$subsystem, Tryptophan Metabolism",),

EC 1.13.11.16 - 2 rhea ids: RHEA:25054 and RHEA:23840
    (rhea_id = , name = "3-carboxyethylcatechol 2,3-dioxygenase", isozymes = [[(4, "A0A1B1EC11"),],], subsystem = subsystem,),

EC 3.7.1.14 - 2 rhea ids: RHEA:34187 and RHEA:34191
(rhea_id = , name = "2-hydroxy-6-oxononatrienedioate hydrolase", isozymes = [[(2, "A0A1B1EC31"),],], subsystem = subsystem,),

already recorded in other subsystems: EC2.6.1.9, EC4.2.1.17, EC1.2.1.10
=#

end # module
