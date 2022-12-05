module Val_Leu_Iso_Degradation

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Valine, Leucine and Isoleucine Degradation"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 32271, name = "3-isopropylmalate dehydrogenase", isozymes = [[(2, "A0A1B1E9B4"),],], subsystem = subsystem,), #1.1.1.85

)

#=
These ECs were present but not included in rxns:

EC 1.1.1.86 - 2 rhea ids: RHEA:22068 and RHEA:13493
    (rhea_id = , name = "Ketol-acid reductoisomerase", isozymes = [[(1, "A0A1B1E8Y0"),],], subsystem = subsystem,),

EC 2.2.1.6 - unclear information for 1 subunit
    (rhea_id = 25249, name = "Acetolactate synthase", isozymes = [[(1, "A0A1B1E9C7"),], [(1, "A0A1B1EG91"),], [(),],], subsystem = subsystem,),
    A0A1B1EA09  SUBUNIT: Dimer of large and small chains. {ECO:0000256|ARBA:ARBA00011744, ECO:0000256|RuleBase:RU368092}

EC 4.2.1.9 - 2 rhea ids: RHEA:24809 and RHEA:27694
    (rhea_id = , name = "Dihydroxy-acid dehydratase", isozymes = [[(2, "A0A1B1EG86"),],], subsystem = subsystem,),

EC 4.2.1.33 - 2 Heterodimers
    (rhea_id = 32287, name = "3-isopropylmalate dehydratase", isozymes = [[(),],], subsystem = subsystem,),
    A0A1B1E9C0  SUBUNIT: Heterodimer of LeuC and LeuD. {ECO:0000256|HAMAP-Rule:MF_01026}
    A0A1B1E9I9  SUBUNIT: Heterodimer of LeuC and LeuD. {ECO:0000256|ARBA:ARBA00011271, ECO:0000256|HAMAP-Rule:MF_01031}

already recorded in other subsystems: EC2.6.1.42, EC4.3.1.19, EC2.3.3.13
=#

end # module
