module Pyruvate_Metabolism

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Pyruvate Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #

rxns = (
    (rhea_id = 21864, name = "Hydroxyacylglutathione hydrolase", isozymes = [[(1, "A0A1B1EEA6"),],], subsystem = subsystem,), #3.1.2.6
    (rhea_id = 19069, name = "Lactoylglutathione lyase", isozymes = [[(1, "A0A1B1EDM9"),],], subsystem = subsystem,), #4.4.1.5
    (rhea_id = 51468, name = "D-lactate dehydrogenase", isozymes = [[(1, "A0A1B1EKM0"),],], subsystem = subsystem,), #1.1.5.12
    (rhea_id = 18157, name = "Pyruvate kinase", isozymes = [[(1, "A0A1B1E9J8"),], [(1, "A0A1B1EDE7"),],], subsystem = "$subsystem, Glycolysis/Gluconeogenesis",), #2.7.1.40
    (rhea_id = 11364, name = "PEP synthase", isozymes = [[(1, "A0A1B1EHZ4"),], [(1, "A0A1B1EK98"),],], subsystem = "$subsystem, Glycolysis/Gluconeogenesis",), #2.7.9.2
    (rhea_id = 28370, name = "Phosphoenolpyruvate carboxylase", isozymes = [[(4, "A0A1B1EFF8"),],], subsystem = subsystem,), #4.1.1.31
    (rhea_id = 15641, name = "Oxaloacetate decarboxylase", isozymes = [[(3, "A0A1B1E9M1"),],], subsystem = subsystem,), #4.1.1.112 - isozyme is also associated with another enzyme (EC4.1.3.17, RHEA:22748, 4-hydroxy-4-methyl-2-oxoglutarate, not present in pyruvate_met but present in vn)?
    (rhea_id = 12653, name = "NAD-dependent malic enzyme", isozymes = [[(4, "A0A1B1ECB3"),],], subsystem = subsystem,), #1.1.1.38
    (rhea_id = 11844, name = "Formate acetyltransferase", isozymes = [[(2, "A0A1B1EB20"),],], subsystem = subsystem,), #2.3.1.54
    (rhea_id = 18181, name = "Malate synthase", isozymes = [[(1, "A0A1B1E9U9"),],], subsystem = subsystem,), #2.3.3.9
    (rhea_id = 11352, name = "Acetate kinase", isozymes = [[(2, "A0A1B1ECS6"),], [(2, "A0A1B1EDJ3"),], [(2, "A0A1B1EJM9"),],], subsystem = subsystem,), #2.7.2.1
    (rhea_id = 19521, name = "Phosphate acetyltransferase", isozymes = [[(6, "A0A1B1EDI1"),], [(1, "A0A1B1EJI5"),],], subsystem = subsystem,), #2.3.1.8
    (rhea_id = 14965, name = "Acylphosphatase", isozymes = [[(1, "A0A1B1EGF9"),],], subsystem = subsystem,), #3.6.1.7
    (rhea_id = 23288, name = "Acetaldehyde dehydrogenase", isozymes = [[(1, "A0A1B1ECS0"),],], subsystem = "$subsystem, Phenylalanine Metabolism",), #1.2.1.10
    (rhea_id = 21524, name = "2-isopropylmalate synthase", isozymes = [[(4, "A0A1B1E9B6"),], [(4, "A0A1B1EE83"),],], subsystem = "$subsystem, Valine, Leucine and Isoleucine Degradation",), #2.3.3.13
    (rhea_id = 12929, name = "Homocitrate synthase", isozymes = [[(1, "A0A1B1EL49"),],], subsystem = "$subsystem, Lysine Biosynthesis",), #2.3.3.14
    (rhea_id = 18617, name = "Phosphoenolpyruvate carboxykinase", isozymes = [[(1, "A0A1B1E8K2"),],], subsystem = "$subsystem, Glycolysis/Gluconeogenesis, Citrate Cycle", ), #4.1.1.49
)

#=
These ECs were present but not included in rxns:

EC 1.2.4.1 - contains a heterodimer
    (rhea_id = 19189, name = "Pyruvate dehydrogenase E1", isozymes = [[(2, "A0A1B1EFH4"),],], subsystem = "$subsystem, Glycolysis/Gluconeogenesis, Citrate Cycle",),
    A0A1B1EFH4  SUBUNIT: Homodimer
    A0A1B1EK05  SUBUNIT: Heterodimer of an alpha and a beta chain. {ECO:0000256|ARBA:ARBA00011870, ECO:0000256|RuleBase:RU361139}
    A0A1B1EK16  SUBUNIT: Heterodimer of an alpha and a beta chain. {ECO:0000256|RuleBase:RU366007}

EC 7.2.4.2 - transporter and includes a heterotrimer
    (rhea_id = 57724, name = "Oxaloacetate decarboxylase", isozymes = [[(),],], subsystem = subsystem,),
    A0A1B1EEW0  SUBUNIT: Heterotrimer of an alpha, a beta and a gamma subunit. {ECO:0000256|ARBA:ARBA00011869, ECO:0000256|HAMAP-Rule:MF_00404}
    A0A1B1EEU7  SUBUNIT: Heterotrimer of an alpha, a beta and a gamma subunit. {ECO:0000256|ARBA:ARBA00011869, ECO:0000256|PIRNR:PIRNR015658}

EC 2.3.1.12 - complicated subunit_structure
    (rhea_id = 17017, name = "Acetyltransferase component of pyruvate dehydrogenase complex", isozymes = [[(),],], subsystem = "$subsystem, Glycolysis/Gluconeogenesis, Citrate Cycle",),
    A0A1B1EES8  SUBUNIT: Forms a 24-polypeptide structural core with octahedral symmetry. {ECO:0000256|ARBA:ARBA00011484, ECO:0000256|RuleBase:RU361137}


These ECs were already recorded in other subsystem files: EC1.8.1.4, EC1.3.5.1, EC4.2.1.2, EC1.1.5.4, EC1.1.1.37, EC 6.2.1.1
=#

end # module
