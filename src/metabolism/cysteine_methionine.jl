module Cys_Met

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Cysteine and Methionine Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 19169, name = "L-serine dehydratase", isozymes = [[(1, "A0A1B1ED13"),],[(1, "A0A1B1EHL2"),],[(1, "A0A1B1EJG8"),],], subsystem = "$subsystem, Glycine, Serine and Threonine Metabolism",), #4.3.1.17
    (rhea_id = 24560, name = "Serine acetyltransferase", isozymes = [[(1, "A0A1B1EFM5"),],], subsystem = subsystem,), #2.3.1.30
    (rhea_id = 14829, name = "Cysteine synthase", isozymes = [[(1, "A0A1B1EAE7"),],], subsystem = subsystem,), #2.5.1.47
    (rhea_id = 12641, name = "D-3-phosphoglycerate dehydrogenase", isozymes = [[(1, "A0A1B1EF12"),],], subsystem = "$subsystem, Glycine, Serine and Threonine Metabolism",), #1.1.1.95
    (rhea_id = 13285, name = "Glutamate--cysteine ligase", isozymes = [[(1, "A0A1B1EEU8"),],], subsystem = subsystem,), #6.3.2.2
    (rhea_id = 13557, name = "Glutathione synthetase", isozymes = [[(1, "A0A1B1EF07"),],], subsystem = subsystem,), #6.3.2.3
    (rhea_id = 22008, name = "Homoserine O-succinyltransferase", isozymes = [[(1, "A0A1B1ECM7"),],], subsystem = subsystem,), #2.3.1.46
    (rhea_id = 24284, name = "Aspartate-semialdehyde dehydrogenase", isozymes = [[(2, "A0A1B1EDS3"),], [(2, "A0A1B1EDU6"),],], subsystem = "$subsystem, Lysine Biosynthesis, Glycine, Serine and Threonine Metabolism",), #1.2.1.11
    (rhea_id = 11172, name = "Methionine synthase", isozymes = [[(1, "A0A1B1EFC9"),],], subsystem = subsystem,), #2.1.1.13
    (rhea_id = 21196, name = "Cobalamin-independent methionine synthase", isozymes = [[(1, "A0A1B1EDA0"),],], subsystem = subsystem,), #2.1.1.14
    (rhea_id = 21080, name = "S-adenosylmethionine synthase", isozymes = [[(1, "A0A1B1EF23"),],], subsystem = subsystem,), #2.5.1.6
    (rhea_id = 17753, name = "S-ribosylhomocysteine lyase", isozymes = [[(2, "A0A1B1EEV7"),],], subsystem = subsystem,), #4.4.1.21
)

#=
These ECs were present but not included in rxns:

EC 2.6.1.52 - 2 rhea ids: RHEA:16573 and RHEA:14329
    (rhea_id = , name = "Phosphoserine aminotransferase", isozymes = [[(2, "A0A1B1ECA5"),],], subsystem = "$subsystem, Glycine, Serine and Threonine Metabolism",),

EC 2.6.1.42 - 3 rhea ids: RHEA:24801, RHEA:18321 and RHEA:24813
    (rhea_id = , name = "Branched-chain-amino-acid aminotransferase", isozymes = [[(1, "A0A1B1EK75"),],], subsystem = "$subsystem, Valine, Leucine and Isoleucine Degradation and Biosynthesis",),

EC 1.1.1.3 - 2 rhea ids: RHEA:15757 and RHEA:15761
    (rhea_id = , name = "Homoserine dehydrogenase", isozymes = [[(4, "A0A1B1E9N6"),], [(4, "A0A1B1EFF7"),],], subsystem = "$subsystem, Lysine Biosynthesis, Glycine, Serine and Threonine Metabolism",), # part of a bifunctional protein

EC 2.7.2.4 - unclear subunit information (unsure which of 2 enzymes the subunit structure refers to)
    (rhea_id = 23776, name = "Aspartokinase", isozymes = [[("A0A1B1E9N6"),], [("A0A1B1EFF7"),], [(1, "A0A1B1EEX2"),], [(1, "A0A1B1EF89"),],], subsystem = "$subsystem, Lysine Biosynthesis, Glycine, Serine and Threonine Metabolism",),

EC 3.2.2.9 - 3 rhea ids: RHEA:29859, RHEA:17805 and RHEA:13617
    (rhea_id = , name = "MTA/SAH nucleosidase", isozymes = [[(1, "A0A1B1E9M0"),],], subsystem = subsystem,),

    
These ECs were already recorded in other subsystem files: EC1.1.1.37
=#

end # module
