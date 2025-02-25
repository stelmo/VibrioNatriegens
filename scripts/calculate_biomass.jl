using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
using SBMLFBCModels
import AbstractFBCModels as A
using RheaReactions, DataFrames, DataFramesMeta, CSV
using FASTX, Statistics

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

molar_masses["Isoleucine"] = 131.173
molar_masses["Leucine"] = 131.1736
molar_masses["Lysine"] = 146.1882
molar_masses["Methionine"] = 149.2124
molar_masses["Phenylalanine"] = 165.19
molar_masses["Threonine"] = 119.1197
molar_masses["Tryptophan"] = 204.2262
molar_masses["Valine"] = 117.1469
molar_masses["Arginine"] = 174.2017
molar_masses["Histidine"] = 155.1552
molar_masses["Alanine"] = 89.0935
molar_masses["Asparagine"] = 132.1184
molar_masses["Aspartate"] = 133.1032 
molar_masses["Cysteine"] = 121.159 
molar_masses["Glutamate"] = 147.1299 
molar_masses["Glutamine"] = 146.1451 
molar_masses["Glycine"] = 75.0669
molar_masses["Proline"] = 115.131
molar_masses["Serine"] = 105.093
molar_masses["Tyrosine"] = 181.1894

molar_masses["glycogen"] = 162.1406 # C6H10O5

molar_masses["fa_c14"] = 228.376
molar_masses["fa_c16"] = 284.484
molar_masses["fa_c18"] = 256.430

# DNA
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

# RNA (assume protein and RNA are directly correlated)

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

# lipids umol/FA
fa_c14_0 = 454
fa_c16_1 = 1628
fa_c16_0 = 1263
fa_c18_1 = 469
fa_c18_0 = 91

# group this into c14 (Tetradecanoic acid), c16 (Hexadecanoic acid), and c18 (Octadecanoic acid) FAs, ignore saturation

fa_c14 = fa_c14_0
fa_c16 = fa_c16_0 + fa_c16_1
fa_c18 = fa_c18_0 + fa_c18_1

# calculate overall biomass components
