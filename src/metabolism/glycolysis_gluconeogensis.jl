module Glycolysis_Gluconeogenesis

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Glycolysis/Gluconeogenesis"

# (rhea_id = , name = , isozymes = [[(),],], subsystem = subsystem,),
 
rxns = (
    (rhea_id = 11816, name = "Glucose-6-phosphate isomerase", isozymes = [ [(1, "A0A1B1EFC4"),], ], subsystem = "$subsystem, Pentose Phosphate Pathway",), #5.3.1.9
    (rhea_id = 11064, name = "Fructose-1,6-bisphosphatase", isozymes = [[(4, "A0A1B1ELE3"),], [(4, "A0A1B1E9A1"),],], subsystem = "$subsystem, Pentose Phosphate Pathway",), #3.1.3.11
    (rhea_id = 16109, name = "ATP-dependent 6-phosphofructokinase", isozymes = [[(4, "A0A1B1EFN6"),],], subsystem = "$subsystem, Pentose Phosphate Pathway",), #2.7.1.11
    (rhea_id = 18585, name = "Triosephosphate isomerase", isozymes = [[(2, "A0A1B1E928"),],], subsystem = subsystem,), #5.3.1.1
    (rhea_id = 14801, name = "Phosphoglycerate kinase", isozymes = [[(1, "A0A1B1EF14"),],], subsystem = subsystem,), #2.7.2.3
    (rhea_id = 14729, name = "Fructose-bisphosphate aldolase", isozymes = [[(1, "A0A1B1EEZ8"),],], subsystem = "$subsystem, Pentose Phosphate Pathway",), #4.1.2.13
    (rhea_id = 17825, name = "Glucokinase", isozymes = [[(1, "A0A1B1EJK5"),],], subsystem = subsystem,), #2.7.1.2
    (rhea_id = 10264, name = "Aldose 1-epimerase", isozymes = [[(1, "A0A1B1ELF1"),],], subsystem = subsystem,), #5.1.3.3
    (rhea_id = 16249, name = "Putative glucose-6-phosphate 1-epimerase", isozymes = [[(1, "A0A1B1EDU3"),],], subsystem = subsystem,), #5.1.3.15
    (rhea_id = 15901, name = "2,3-bisphosphoglycerate-independent phosphoglycerate mutase", isozymes = [[(1, "A0A1B1EFI9"),],], subsystem = "$subsystem, Glycine, Serine and Threonine Metabolism",), #5.4.2.12
    (rhea_id = 15045, name = "Dihydrolipoyl dehydrogenase", isozymes = [[(1, "A0A1B1EET9"),],], subsystem = "$subsystem, Citrate Cycle, Pyruvate Metabolism, Glycine, Serine and Threonine Metabolism, Valine, Leucine and Isoleucine Biosynthesis",), # 1.8.1.4
    (rhea_id = 23176, name = "AcCoA synthetase", isozymes = [[(1, "A0A1B1EFQ7"),],], subsystem = "$subsystem, Pyruvate Metabolism",), #6.2.1.1
)


#=
These ECs were present but not included in rxns:
EC 2.7.1.199 - PTS transporter
    (rhea_id = , name = , isozymes = [[(),],], subsystem = subsystem,),

EC 4.2.1.11 - unclear subunit information
    (rhea_id = , name = , isozymes = [[(),],], subsystem = subsystem,),


These ECs were already recorded in other subsystem files: EC 2.3.1.12, EC 5.4.2.2, EC 4.1.1.49, EC 1.2.4.1, EC 2.7.1.40, EC 2.7.9.2
=#

end # module
