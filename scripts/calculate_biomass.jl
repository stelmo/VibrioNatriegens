using JSON
using DataFrames, DataFramesMeta, CSV
using FASTX, Statistics

biomass = Dict()

chebi_lookup = Dict(
    "glycogen" => "glycogen",

    "GTP" => "CHEBI:37565", # GTP
    "CTP" => "CHEBI:37563", # CTP
    "UTP" => "CHEBI:46398", # UTP
    "ATP" => "CHEBI:30616", # atp

    "dGTP" => "CHEBI:61429", # dGTP
    "dATP" => "CHEBI:61404", # dATP
    "dCTP" => "CHEBI:61481", # dCTP
    "dTTP" => "CHEBI:37568", # dTTP

    "K" => "CHEBI:32551", # lysine
    "M" => "CHEBI:57844", # methionine
    "C" => "CHEBI:35235", # cysteine
    "A" => "CHEBI:57972", # alanine
    "D" => "CHEBI:29991", # aspartate
    "N" => "CHEBI:58048", # asparagine
    "E" => "CHEBI:29985", # glutamate
    "Q" => "CHEBI:58359", # glutamine
    "T" => "CHEBI:57926", # threonine
    "S" => "CHEBI:33384", # serine
    "G" => "CHEBI:57305", # glycine
    "I" => "CHEBI:58045", # isoleucine
    "V" => "CHEBI:57762", # valine
    "L" => "CHEBI:57427", # leucine
    "R" => "CHEBI:32682", # arginine
    "H" => "CHEBI:57595", # histidine
    "Y" => "CHEBI:58315", # tyrosine
    "F" => "CHEBI:58095", # phenylalanine
    "W" => "CHEBI:57912", # tryptophan
    "P" => "CHEBI:60039", # proline

    "tetra" => "CHEBI:30807", # tetradecanoic acid
    "hexa" => "CHEBI:7896", # hexadecanoic acid
    "octa" => "CHEBI:25629", # octadecanoic acid   
)

# Long 2017 data g/g
protein = 0.465
rna = 0.286
lipid = 0.075
glycogen = 0.034

# e coli g/g
dna = 0.01 # 87% up to here

# collate all molar masses
molar_masses = Dict() # g/mol

molar_masses["dATP"] = 491.181
molar_masses["dGTP"] = 507.181
molar_masses["dCTP"] = 467.156
molar_masses["dTTP"] = 482.168

molar_masses["ATP"] = 507.18
molar_masses["GTP"] = 523.180
molar_masses["CTP"] = 483.156
molar_masses["UTP"] = 484.139

molar_masses["I"] = 131.173
molar_masses["L"] = 131.1736
molar_masses["K"] = 146.1882
molar_masses["M"] = 149.2124
molar_masses["F"] = 165.19
molar_masses["T"] = 119.1197
molar_masses["W"] = 204.2262
molar_masses["V"] = 117.1469
molar_masses["R"] = 174.2017
molar_masses["H"] = 155.1552
molar_masses["A"] = 89.0935
molar_masses["N"] = 132.1184
molar_masses["D"] = 133.1032 
molar_masses["C"] = 121.159 
molar_masses["E"] = 147.1299 
molar_masses["Q"] = 146.1451 
molar_masses["G"] = 75.0669
molar_masses["P"] = 115.131
molar_masses["S"] = 105.093
molar_masses["Y"] = 181.1894

molar_masses["glycogen"] = 162.1406 # C6H10O5

molar_masses["tetra"] = 228.376
molar_masses["hexa"] = 284.484
molar_masses["octa"] = 256.430

# DNA
dna_lookup = Dict(
   'A' => "dATP",
   'T' => "dTTP",  
   'G' => "dGTP",  
   'C' => "dCTP",  
)
dna_bases = Dict()
FASTAReader(open(joinpath("data", "genome", "genome.fasta"))) do reader
    for record in reader
    for s in sequence(record)
        dna_bases[s] = get(dna_bases, s, 0) + 1
    end
    end
end
dna_total = sum(values(dna_bases))
for (k, v) in dna_bases
    dna_bases[k] = v / dna_total
end
for (k, v) in dna_bases
    biomass[chebi_lookup[dna_lookup[k]]] = -v * dna / molar_masses[dna_lookup[k]] * 1000
end

# proteome
df = DataFrame(CSV.File("quantitative_proteome.csv")) # use glucose as the baseline
@subset!(df, :Condition .== "Glucose")
df = @combine(groupby(df, :Protein), :MeanMoleFraction = mean(:MoleFraction))

plu = Dict(df.Protein .=> df.MeanMoleFraction)

protein_bases = Dict()
FASTAReader(open(joinpath("data", "genome", "proteome.fasta"))) do reader
    for record in reader
    desc = description(record)
    if occursin("protein_id=", desc)
        pid = string(first(split(string(last(split(desc, "protein_id="))),"]"))) 
        for s in sequence(record)
        protein_bases[s] = get(protein_bases, s, 0) + 1 * get(plu, pid, 0.0)
        end
    end
    end
end

protein_total = sum(values(protein_bases))
for (k, v) in protein_bases
    protein_bases[k] = v / protein_total
end
for (k, v) in protein_bases
    biomass[chebi_lookup[string(k)]] = -v * protein / molar_masses[string(k)] * 1000
end

# RNA (assume protein and RNA are directly correlated)
rna_lookup = Dict(
    'A' => "ATP",
    'U' => "UTP",  
    'G' => "GTP",  
    'C' => "CTP",
)
rna_bases = Dict()
FASTAReader(open(joinpath("data", "genome", "transcriptome.fasta"))) do reader
    for record in reader
    desc = description(record)
    if occursin("protein_id=", desc)
        pid = string(first(split(string(last(split(desc, "protein_id="))),"]"))) 
        for s in sequence(record)
        rna_bases[s] = get(rna_bases, s, 0) + 1 * get(plu, pid, 0.0)
        end
    end
    end
end

rna_total = sum(values(rna_bases))
for (k, v) in rna_bases
    rna_bases[k] = v / rna_total
end
rna_bases['U'] = rna_bases['T']
delete!(rna_bases, 'T')
for (k, v) in rna_bases
    biomass[chebi_lookup[rna_lookup[k]]] = -v * protein / molar_masses[rna_lookup[k]] * 1000
end

# lipids umol/FA
fa_c14_0 = 454
fa_c16_1 = 1628
fa_c16_0 = 1263
fa_c18_1 = 469
fa_c18_0 = 91

# group this into c14 (Tetradecanoic acid), c16 (Hexadecanoic acid), and c18 (Octadecanoic acid) FAs, ignore saturation
fa = Dict(
    "tetra" => fa_c14_0,
    "hexa" =>  fa_c16_0 + fa_c16_1,
    "octa" => fa_c18_0 + fa_c18_1,
)

for (k, v) in fa
    biomass[chebi_lookup[k]] = -v * lipid * 0.001 # to mmol/gDW
end

# glycogen
biomass["glycogen"] = -glycogen * 1000 / molar_masses["glycogen"]


# required atp for growth
atp_req = 75.0

biomass["CHEBI:30616"] = biomass["CHEBI:30616"] - atp_req # atp
biomass["CHEBI:15377"] = -atp_req # water
biomass["CHEBI:43474"] = atp_req # pi
biomass["CHEBI:456216"] = atp_req # adp
biomass["CHEBI:15378"] = atp_req # h+

biomass

open(joinpath("data", "model", "biomass.json"), "w") do io
    JSON.print(io, biomass)
end
