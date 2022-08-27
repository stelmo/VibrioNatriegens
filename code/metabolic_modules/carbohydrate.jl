using Revise

includet("code//reconstruction.jl")
import .Reconstruction as rc

using COBREXA, Tulip, Clarabel 
using DataFrames, DataFramesMeta, CSV, Cleaner
using KEGGModuleParser
import BiGGReactions as br
import RheaReactions as rr
import MetaNetXReactions as mnr

vibrio_df = DataFrame(CSV.File(joinpath("data", "uniprot", "v_natriegens.tsv"))) |> CleanTable |> polish_names! |> DataFrame
names(vibrio_df)
vibrio_df[!, "e_c_number"]

model = StandardModel("Vibrio_Natriegens")

#: ED: https://www.genome.jp/pathway/map00030+M00008
# ecs = rc.get_ecs("M00008")
# df = @subset vibrio_df @byrow begin 
#     !ismissing(:e_c_number) && any(in.(:e_c_number, Ref(ecs)))
# end
# select!(df, [:entry, :protein_names, :e_c_number, :rhea_i_d, :subunit_structure])
# df
# mnr.get_reaction_from_rhea(17089).crossreferences

rc.add_reaction_from_rhea!(
    model, 
    15841; 
    name="G6PDH2c", 
    grrs=[["A0A1B1ECJ9"]],
    subsystem="Entner-Doudoroff pathway",
)

rc.add_reaction_from_rhea!(
    model, 
    17277; 
    name="EDD", 
    grrs=[["A0A1B1E888"]],
    subsystem="Entner-Doudoroff pathway",
)

rc.add_reaction_from_rhea!(
    model, 
    12556; 
    name="PGL", 
    grrs=[["A0A1B1ECI2"]],
    subsystem="Entner-Doudoroff pathway",
)

rc.add_reaction_from_rhea!(
    model, 
    17089; 
    name = "KDPGALDOL-RXN",
    grrs=[["A0A1B1E8D0"], ["A0A1B1EGV6"], ["A0A1B1EHS6",]],
    subsystem="Entner-Doudoroff pathway",
)

model

# EMP 
ecs = rc.get_ecs("M00001")
df = @subset vibrio_df @byrow begin 
    !ismissing(:e_c_number) && any(in.(:e_c_number, Ref(ecs)))
end
select!(df, [:entry, :protein_names, :e_c_number, :rhea_i_d, :subunit_structure])
df
mnr.get_reaction_from_rhea(17089).crossreferences