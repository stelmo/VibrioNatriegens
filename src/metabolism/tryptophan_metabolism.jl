module Try_Metabolism

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Tryptophan Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 19553, name = "Tryptophanase", isozymes = [[(4, "A0A1B1EI58"),],], subsystem = subsystem,), #4.1.99.1
    (rhea_id = 52848, name = "tRNA U34 carboxymethyltransferase", isozymes = [[(4, "A0A1B1EBE5"),],], subsystem = subsystem,), #2.5.1.-

)

#=
These ECs were present but not included in rxns:

EC 2.6.1.- missing rhea id
    (rhea_id = , name = "Aminotransferase", isozymes = [[(2, "A0A1B1ED36"),], [(2, "A0A1B1ELW5"),],], subsystem = subsystem,),

EC 1.2.1.- missing rhea id
    (rhea_id = , name = "Glyceraldehyde-3-phosphate dehydrogenase", isozymes = [[(1, "A0A1B1EDV5"),], [(1, "A0A1B1EFY0"),], [(1, "A0A1B1EH68"),], [(1, "A9X8D5"),], [(1, "S5FTK5"),],], subsystem = subsystem,),

EC 1.3.1.- 4 rhea ids: RHEA:54452, RHEA:23624, RHEA:53380 and RHEA:53376
    (rhea_id = , name = "tRNA-dihydrouridine synthase", isozymes = [[(1, "A0A1B1EDT4"),], [(1, "A0A1B1EGN8"),],], subsystem = subsystem,),

EC 6.3.2.- missing rhea id
    (rhea_id = , name = "Alpha-L-glutamate ligase", isozymes = [[(1, "A0A1B1ED93"),], [(1, "A0A1B1EHM5"),],], subsystem = subsystem,),


already recorded in other subsystems: EC2.1.1.-, EC1.11.1.21, EC1.1.1.35, EC4.2.1.7 and EC1.2.4.2
=#

end # module
