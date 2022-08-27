function carbohydrate!(model)
#=
This script builds the carbohydrate metabolism. The import statements were used
to identify reactions to build the model, they should be uncommented only when
troubleshooting.
=#
# using Revise
# includet("code//reconstruction.jl")
# import .Reconstruction as rc
# using COBREXA, Tulip, Clarabel 
# using DataFrames, DataFramesMeta, CSV, Cleaner
# using KEGGModuleParser
# import BiGGReactions as br
# import RheaReactions as rr
# import MetaNetXReactions as mnr

# vibrio_df = DataFrame(CSV.File(joinpath("data", "uniprot", "v_natriegens.tsv"))) |> CleanTable |> polish_names! |> DataFrame
# names(vibrio_df)
# vibrio_df[!, "e_c_number"]

# model = StandardModel("Vibrio_Natriegens")

#=
#: Entner-Doudoroff pathway 
https://www.genome.jp/pathway/map00030+M00008

Uncomment the dataframe below only for troubleshooting. Notes stores the 
grr stoichiometry if not monomers or monomeric complexes.
=#
# ecs = get_ecs("M00008")
# df = @subset vibrio_df @byrow begin 
#     !ismissing(:e_c_number) && any(in.(:e_c_number, Ref(ecs)))
# end
# select!(df, [:entry, :protein_names, :e_c_number, :rhea_i_d, :subunit_structure])
# df
# mnr.get_reaction_from_rhea(17089).crossreferences

ed_pathway = [
    (15841, "G6PDH2c", [["A0A1B1ECJ9"]], Notes()),
    (17277, "EDD", [["A0A1B1E888"]], Notes()),
    (12556, "PGL", [["A0A1B1ECI2"]], Notes()),
    (17089, "KHG/KDPG aldolase", [["A0A1B1E8D0"], ["A0A1B1EGV6"], ["A0A1B1EHS6",]], note_complex("[[3],[3],[3]")), # gluconate induced, homotrimeric in e coli
]

for (rhea_id, bigg_id, grrs, notes) in ed_pathway
    add_reaction_from_rhea!(
        model, 
        rhea_id; 
        name=bigg_id, 
        grrs,
        subsystem="Entner-Doudoroff pathway",
        notes,
    )
end

#=
Glycolysis (Embden-Meyerhof pathway)
https://www.genome.jp/pathway/map00010+M00001

Uncomment the dataframe below only for troubleshooting. Notes stores the 
grr stoichiometry if not monomers or monomeric complexes.
=# 
# ecs = get_ecs("M00001")
# df = @subset vibrio_df @byrow begin 
#     !ismissing(:e_c_number) && any(in.(:e_c_number, Ref(ecs)))
# end
# select!(df, [:entry, :protein_names, :e_c_number, :rhea_i_d, :subunit_structure])
# df
# mnr.get_reaction_from_rhea(17825).crossreferences

emp_pathway = [
    (16109, "PFK", [["A0A1B1EFN6"]], note_complex("[[4]]")),
    (18585, "TPI", [["A0A1B1E928"]], note_complex("[[2]]")),
    (14801, "PGK", [["A0A1B1EF14"]], Notes()),
    (11816, "PGI", [["A0A1B1EFC4"]], Notes()),
    (15901, "PGM", [["A0A1B1EFI9"]], Notes()),
    (10164, "ENO", [["A0A1B1EFN0"]], note_complex("[[2]]")),
    (18157, "PYK", [["A0A1B1E9J8"], ["A0A1B1EDE7"]], Notes()),
    (14729, "Fructose-bisphosphate aldolase", [["A0A1B1EEZ8"]], Notes()),
    (17053, "Transaldolase", [["A0A1B1EGX7"]], Notes()),
    (17825, "HEX1", [["A0A1B1EJK5"]], Notes()),
]

for (rhea_id, bigg_id, grrs) in emp_pathway
    add_reaction_from_rhea!(
        model, 
        rhea_id; 
        name=bigg_id, 
        grrs,
        subsystem="Embden-Meyerhof pathway",
    )
end

end # build function