module Purine_Metabolism 

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Purine Metabolism"

rxns = (
    (rhea_id = 17453, name = "Phosphoribosylamine--glycine ligase", isozymes = [[(1, "A0A1B1EFN9" ),],], subsystem = subsystem,),#6.3.4.13
    (rhea_id = 15053, name = "Phosphoribosylglycinamide formyltransferase", isozymes = [[(1, "A0A1B1EE84" ),],], subsystem = subsystem,),#2.1.2.2
    (rhea_id = 24829, name = "Formate-dependent phosphoribosylglycinamide formyltransferase", isozymes = [[(2, "A0A1B1EBR7" ),],], subsystem = subsystem,),#6.3.1.21
    (rhea_id = 17129, name = "Phosphoribosylformylglycinamidine synthase", isozymes = [[(1, "A0A1B1EA44" ),],], subsystem = subsystem,),#6.3.5.3
    (rhea_id = 23032, name = "Phosphoribosylformylglycinamidine cyclo-ligase", isozymes = [[(1, "A0A1B1EE96" ),],], subsystem = subsystem,),#6.3.3.1
    (rhea_id = 19317, name = "N5-carboxyaminoimidazole ribonucleotide synthase", isozymes = [[(2, "A0A1B1EG27" ),],], subsystem = subsystem,),#6.3.4.18
    (rhea_id = 13193, name = "N5-carboxyaminoimidazole ribonucleotide mutase", isozymes = [[(1, "A0A1B1EG13" ),],], subsystem = subsystem,),#5.4.99.18
    (rhea_id = 22628, name = "Phosphoribosylaminoimidazole-succinocarboxamide synthase", isozymes = [[(1, "A0A1B1EBP3" ),],], subsystem = subsystem,),#6.3.2.6
    (rhea_id = 22192, name = "Bifunctional purine biosynthesis protein PurH", isozymes = [[(1, "A0A1B1EG16" ),],], subsystem = subsystem,),#2.1.2.3 cen
    (rhea_id = 10412, name = "ADP-ribose pyrophosphatase", isozymes = [[(1, "A0A1B1EA47" ),],], subsystem = subsystem,),#3.6.1.13
    (rhea_id = 18445, name = "Bifunctional purine biosynthesis protein PurH", isozymes = [[(1, "A0A1B1EG16" ),],], subsystem = subsystem,),#3.5.4.10 ?
    (rhea_id = 24252, name = "Bis(5'-nucleosyl)-tetraphosphatase", isozymes = [[(1, "A0A1B1E9A2" ),],], subsystem = subsystem,),#3.6.1.41
    (rhea_id = 24152, name = "Adenylyl-sulfate kinase", isozymes = [[(1, "A0A1B1E986" ),],], subsystem = subsystem,),#2.7.1.25
    (rhea_id = 21528, name = "Exopolyphosphatase", isozymes = [[(1, "A0A1B1E9U1" ),],], subsystem = subsystem,),#3.6.1.11
    (rhea_id = 13073, name = "Guanosine-5'-triphosphate,3'-diphosphate pyrophosphatase", isozymes = [[(1, "A0A1B1EFZ9" ),],], subsystem = subsystem,),#3.6.1.40
    (rhea_id = 14253, name = "Guanosine-3',5'-bis(diphosphate) 3'-diphosphatase", isozymes = [[(1, "A0A1B1E8Q7" ),],], subsystem = subsystem,),#3.1.7.2
    (rhea_id = 25277, name = "3',5'-cyclic adenosine monophosphate phosphodiesterase CpdA", isozymes = [[(1, "A0A1B1E9G4" ),],], subsystem = subsystem,),#3.1.4.53
    (rhea_id = 20780, name = "Guanylate kinase", isozymes = [[(1, "A0A1B1E9A8" ),],], subsystem = subsystem,),#2.7.4.8
    (rhea_id = 17973, name = "Hypoxanthine phosphoribosyltransferase", isozymes = [[(1, "A0A1B1EER0" ),],], subsystem = subsystem,),#2.4.2.8
    (rhea_id = 11680, name = "GMP synthase [glutamine-hydrolyzing]", isozymes = [[(2, "A0A1B1E9Z6" ),],], subsystem = subsystem,),#6.3.5.2
    (rhea_id = 11708, name = "Inosine-5'-monophosphate dehydrogenase", isozymes = [[(4, "A0A1B1E9Y2" ),],], subsystem = subsystem,),#1.1.1.205
    (rhea_id = 14665, name = "tRNA-specific adenosine deaminase", isozymes = [[(1, "A0A1B1ELM7" ),],], subsystem = subsystem,),#3.5.4.3
    (rhea_id = 23688, name = "Adenine deaminase", isozymes = [[(1, "A0A1B1EH73" ),],], subsystem = subsystem,),#3.5.4.2
    (rhea_id = 20129, name = "AMP nucleosidase", isozymes = [[(1, "A0A1B1EA82" ),],], subsystem = subsystem,),#3.2.2.4
    (rhea_id = 16609, name = "Adenine phosphoribosyltransferase", isozymes = [[(2, "A0A1B1EE60" ),],], subsystem = subsystem,),#2.4.2.7
    (rhea_id = 12973, name = "Adenylate kinase", isozymes = [[(1, "A0A1B1EAG5" ),],], subsystem = subsystem,),#2.7.4.3
    (rhea_id = 23736, name = "5-hydroxyisourate hydrolase", isozymes = [[(4, "A0A1B1EHA3" ),],], subsystem = subsystem,),#3.5.2.17
)


#=
These ECs were present but not included in rxns:
(rhea_id = 27658/18793, name = "Phosphopentomutase", isozymes = [[(1, "A0A1B1EEL4"),],], subsystem = subsystem,), #5.4.2.7
(rhea_id = ?, name = "dITP/XTP pyrophosphatase", isozymes = [[(2, "A0A1B1EF11" ),],], subsystem = subsystem,),#3.6.1.66 missing reha
(rhea_id = 18133, name = "Sulfate adenylyltransferase subunit 1", isozymes = [[(2, "A0A1B1E975" ),], [(2, "A0A1B1E9E5"),],], subsystem = subsystem,),#2.7.7.4 heterodimer
(rhea_id = 28330/28406 , name = "Inosine/xanthosine triphosphatase", isozymes = [[(2, "A0A1B1E9Z8" ),],], subsystem = subsystem,),#3.6.1.73
# strange information for 2.4.2.1
(rhea_id = 25424/17973/10800 , name = "Xanthine-guanine phosphoribosyltransferase", isozymes = [[(4, "A0A1B1EA25" ),],], subsystem = subsystem,),#2.4.2.22
(rhea_id = 27710/21140, name = "Guanosine-inosine kinase", isozymes = [[(1, "A0A1B1EK80" ),], [(1, "A0A1B1EBD4"),],], subsystem = subsystem,),#2.7.1.73
(rhea_id = 28190/24408, name = "Adenosine deaminase", isozymes = [[(1, "A0A1B1E8H2" ),],], subsystem = subsystem,),#3.5.4.4
=#

end # module
