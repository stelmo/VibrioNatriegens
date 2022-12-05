module Val_Leu_Iso_Synthesis

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Valine, Leucine and Isoleucine Biosynthesis"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 14477, name = "Acyl-coenzyme A dehydrogenase", isozymes = [[(1, "A0A1B1EBK8"),], [(1, "A0A1B1EE66"),],], subsystem = subsystem,), #1.3.8.7
    (rhea_id = 24404, name = "Hydroxymethylglutaryl-CoA lyase", isozymes = [[(1, "A0A1B1EJY5"),],], subsystem = subsystem,), #4.1.3.4
    (rhea_id = 17681, name = "3-hydroxyisobutyrate dehydrogenase", isozymes = [[(1, "A0A1B1EJU5"),], [(1, "A0A1B1ELB3"),],], subsystem = subsystem,), #1.1.1.31
)

#=
These ECs were present but not included in rxns:

EC 2.3.1.16 - Heterotetramer
    (rhea_id = 21564, name = "3-ketoacyl-CoA thiolase", isozymes = [[(),],], subsystem = subsystem,),
    A0A1B1E8F4   SUBUNIT: Heterotetramer of two alpha chains (FadB) and two beta chains (FadA). {ECO:0000256|HAMAP-Rule:MF_01620}

already recorded in other subsystems: EC2.6.1.42, EC1.8.1.4, EC4.2.1.7 and EC1.1.1.35
=#

end # module
