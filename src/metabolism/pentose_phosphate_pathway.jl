module Pentose_Phosphate_Pathway

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Pentose Phosphate Pathway"

# (rhea_id = , name = , isozymes = [[(),],], subsystem = subsystem,),
 
rxns = (
    (rhea_id = 19433, name = "Gluconokinase", isozymes = [[(1, "A0A1B1E8C0"),], [(1, "A0A1B1EL04"),],], subsystem = subsystem,), #EC 2.7.1.12
    (rhea_id = 17277, name = "MurNAc-6-P etherase", isozymes = [[(1, "A0A1B1E888"),],], subsystem = subsystem,), #4.2.1.12
    (rhea_id = 17089, name = "2-dehydro-3-deoxy-phosphogluconate aldolase", isozymes = [[(1, "A0A1B1E8D0"),], [(1, "A0A1B1EGV6"),], [(1, "A0A1B1EHS6"),],], subsystem = subsystem,), #EC 4.1.2.14
    (rhea_id = 12556, name = "6-phosphogluconolactonase", isozymes = [[(1, "A0A1B1ECI2"),],], subsystem = subsystem,), #3.1.1.31
    (rhea_id = 15841, name = "Glucose-6-phosphate 1-dehydrogenase", isozymes = [[(1, "A0A1B1ECJ9"),],], subsystem = subsystem,), #1.1.1.49
    (rhea_id = 10116, name = "6-phosphogluconate dehydrogenase", isozymes = [[(2, "A0A1B1ECK0"),],], subsystem = subsystem,), #1.1.1.44
    (rhea_id = 13677, name = "Ribulose-phosphate 3-epimerase", isozymes = [[(1, "A0A1B1EFE9"),],], subsystem = subsystem,), #5.1.3.1
    (rhea_id = 10508, name = "Transketolase", isozymes = [[(2, "A0A1B1EEZ3"),], [(2, "A0A1B1EGZ9"),],], subsystem = subsystem,), #2.2.1.1
    (rhea_id = 17053, name = "Transaldolase", isozymes = [[(1, "A0A1B1EGX7"),],], subsystem = subsystem,), #2.2.1.2
    (rhea_id = 14657, name = "Ribose-5-phosphate isomerase A", isozymes = [[(2, "A0A1B1EF17"),],], subsystem = subsystem,), #5.3.1.6
    (rhea_id = 15609, name = "Ribose-phosphate pyrophosphokinase", isozymes = [[(6, "A0A1B1EA88"),],], subsystem = "$subsystem, Purine Metabolism"), #2.7.6.1
    (rhea_id = 13697, name = "Ribokinase", isozymes = [[(2, "A0A1B1EAH9"),], [(2, "A0A1B1EJU7"),], [(2, "A0A1B1ELP3"),],], subsystem = subsystem,), #EC 2.7.1.15
    (rhea_id = 23536, name = "Phosphoglucomutase", isozymes = [[(1, "A0A1B1EB72"),], [(1, "A0A1B1EER5"),],], subsystem = "$subsystem, Glycolysis/Gluconeogenesis, Purine Metabolism" ), #EC 5.4.2.2
    (rhea_id = 12821, name = "Deoxyribose-phosphate aldolase", isozymes = [[(1, "A0A1B1EEJ0"),],], subsystem = subsystem,), #4.1.2.4
)

#=
These ECs were present but not included in rxns:
EC 5.4.2.7 - 2 rhea ids: 27658 and 18793
    (rhea_id = , name = "Phosphopentomutase", isozymes = [[(1, "A0A1B1EEL4"),],], subsystem = "$subsystem, Purine Metabolism",),

These ECs were already recorded in other subsystem files: EC5.3.1.9, EC3.1.3.11, EC2.7.1.11, EC4.1.2.13
=#

end # module
