using JSON
using DataFrames, DataFramesMeta, CSV
using FASTX, Statistics

biomass = Dict()

chebi_lookup = Dict(
    "glycogen" => "glycogen",
    "GTP" => "37565", # GTP
    "CTP" => "37563", # CTP
    "UTP" => "46398", # UTP
    "ATP" => "30616", # atp
    "dGTP" => "61429", # dGTP
    "dATP" => "61404", # dATP
    "dCTP" => "61481", # dCTP
    "dTTP" => "37568", # dTTP
    "K" => "32551", # lysine
    "M" => "57844", # methionine
    "C" => "35235", # cysteine
    "A" => "57972", # alanine
    "D" => "29991", # aspartate
    "N" => "58048", # asparagine
    "E" => "29985", # glutamate
    "Q" => "58359", # glutamine
    "T" => "57926", # threonine
    "S" => "33384", # serine
    "G" => "57305", # glycine
    "I" => "58045", # isoleucine
    "V" => "57762", # valine
    "L" => "57427", # leucine
    "R" => "32682", # arginine
    "H" => "57595", # histidine
    "Y" => "58315", # tyrosine
    "F" => "58095", # phenylalanine
    "W" => "57912", # tryptophan
    "P" => "60039", # proline
    "tetra" => "30807", # tetradecanoic acid
    "hexa" => "7896", # hexadecanoic acid
    "octa" => "25629", # octadecanoic acid   
)

# Long 2017 data g/g
protein = 0.458
rna = 0.286
lipid = 0.075
glycogen = 0.034

# e coli g/g
peptidoglycan = 0.025
dna = 0.031
lipopolysaccharide = 0.034

total = protein + rna + lipid + glycogen + dna + peptidoglycan + lipopolysaccharide

soluble_pool = 1.0 - total # overestimate

# collate all molar masses
molar_masses = Dict() # g/mol
begin
    molar_masses["61404"] = 487.1499 # dATP
    molar_masses["61429"] = 503.1493 # dGTP
    molar_masses["61481"] = 463.1252 # dCTP
    molar_masses["37568"] = 478.1365 # dTTP
    
    molar_masses["30616"] = 503.14946 # ATP
    molar_masses["37565"] = 519.14886 # GTP
    molar_masses["37563"] = 479.12468 # CTP
    molar_masses["46398"] = 480.1094 # UTP

    molar_masses["32551"] = 147.19558 # # lysine
    molar_masses["58045"] = 131.1729 # isoleucine
    molar_masses["57427"] = 131.1729 # leucine
    molar_masses["57844"] = 149.2124 # methionine
    molar_masses["58095"] = 165.1891 # phenylalanine
    molar_masses["57926"] = 119.1197  # threonine
    molar_masses["57912"] = 204.2262 # tryptophan
    molar_masses["57762"] = 117.1469 # valine
    molar_masses["32682"] = 175.20906 # arginine
    molar_masses["57595"] = 155.1552 # histidine
    molar_masses["57972"] = 89.0935 # alanine
    molar_masses["58048"] = 132.1184 # asparagine
    molar_masses["29991"] = 132.09478 # aspartate
    molar_masses["35235"] = 121.159 # cysteine
    molar_masses["29985"] = 147.1299 # glutamate
    molar_masses["58359"] = 146.1451 # glutamine
    molar_masses["57305"] = 75.0669 # glycine
    molar_masses["60039"] = 115.131 # proline
    molar_masses["33384"] = 105.093 # serine
    molar_masses["58315"] = 181.1894 # tyrosine

    molar_masses["glycogen"] = 162.1406 # C6H10O5
  
    molar_masses["30807"] = 227.364 # tetradecanoic acid
    molar_masses["7896"] = 255.4161 # hexadecanoic acid
    molar_masses["25629"] = 283.47 # octadecanoic acid

    molar_masses["peptidoglycan"] = 1916.20990

    molar_masses["kdo_lps"] = 2232.67080
    
    # soluble pool
    molar_masses["60530"] = 836.838 # Fe(II)-heme o
    molar_masses["57692"] = 782.5259 # FAD
    molar_masses["57705"] = 605.3378 # UDP-N-acetyl-alpha-D-glucosamine
    molar_masses["57540"] = 662.4172 # NAD(+)
    molar_masses["58885"] = 564.2859 # UDP-alpha-D-glucose
    molar_masses["57287"] = 763.502 # CoA
    molar_masses["57925"] = 306.31 # glutathione
    molar_masses["57945"] = 663.4251 # NADH
    molar_masses["58223"] = 401.1374 # UDP
    molar_masses["29985"] = 146.12136 # L-glutamate
    molar_masses["32966"] = 336.08392 # beta-D-fructose 1,6-bisphosphate
    molar_masses["30616"] = 503.14946 # ATP
    molar_masses["57783"] = 741.3891 # NADPH
    molar_masses["57986"] = 375.356 # riboflavin
    molar_masses["597326"] = 245.126 # pyridoxal 5'-phosphate
    molar_masses["62501"] = 439.3816 # folate
    molar_masses["58297"] = 610.615 # glutathione disulfide
    molar_masses["58210"] = 453.325 # FMN
    molar_masses["58349"] = 740.3812 # NADP(+)
end

for (k, _) in model.reactions["biomass"].stoichiometry
    haskey(molar_masses, k) || continue
    println(k, ": ", molar_masses[k] -  parse(Float64, first(model.metabolites[k].annotations["molarmass"])))
end

# DNA
dna_lookup = Dict('A' => "dATP", 'T' => "dTTP", 'G' => "dGTP", 'C' => "dCTP")
dna_bases = Dict()
FASTAReader(open(joinpath("data", "annotations", "genome", "genome.fasta"))) do reader
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
    ch = chebi_lookup[dna_lookup[k]]
    biomass[ch] = -v * dna / molar_masses[ch] * 1000
end

# proteome
df = DataFrame(CSV.File("quantitative_proteome.csv")) # use glucose as the baseline
@subset!(df, :Condition .== "Glucose")
df = @combine(groupby(df, :Protein), :MeanMoleFraction = mean(:MoleFraction))

plu = Dict(df.Protein .=> df.MeanMoleFraction)

protein_bases = Dict()
FASTAReader(open(joinpath("data", "annotations", "genome", "proteome.fasta"))) do reader
    for record in reader
        desc = description(record)
        if occursin("protein_id=", desc)
            pid = string(first(split(string(last(split(desc, "protein_id="))), "]")))
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
    ch = chebi_lookup[string(k)]
    biomass[ch] = -v * protein / molar_masses[ch] * 1000
end

# RNA (assume protein and RNA are directly correlated)
rna_lookup = Dict('A' => "ATP", 'U' => "UTP", 'G' => "GTP", 'C' => "CTP")
rna_bases = Dict()
FASTAReader(open(joinpath("data", "annotations", "genome", "transcriptome.fasta"))) do reader
    for record in reader
        desc = description(record)
        if occursin("protein_id=", desc)
            pid = string(first(split(string(last(split(desc, "protein_id="))), "]")))
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
    ch = chebi_lookup[rna_lookup[k]]
    biomass[ch] = -v * rna / molar_masses[ch] * 1000
end

# lipids umol/FA
fa_c14_0 = 454
fa_c16_1 = 1628
fa_c16_0 = 1263
fa_c18_1 = 469
fa_c18_0 = 91

# group this into c14 (Tetradecanoic acid), c16 (Hexadecanoic acid), and c18 (Octadecanoic acid) FAs, ignore saturation
fa = Dict("tetra" => fa_c14_0, "hexa" => fa_c16_0 + fa_c16_1, "octa" => fa_c18_0 + fa_c18_1)

for (k, v) in fa
    biomass[chebi_lookup[k]] = -v * lipid * 0.001 # to mmol/gDW
end

# glycogen
biomass["glycogen"] = -glycogen * 1000 / molar_masses["glycogen"]

biomass["61388"] = -peptidoglycan * 1000 / molar_masses["peptidoglycan"]

biomass["58540"] = -lipopolysaccharide * 1000 / molar_masses["kdo_lps"]

# add vitamins here at a very small fraction
solubles = Dict( # molar fractions
    "29985" => 10.0, # glutamate
    "57925" => 2.0, # glutathione
    "32966" => 2.0, # fructose bisphosphate
    "30616" => 1.0, # atp
    "57705" => 1.0, # UDP-N-acetyl-alpha-D-glucosamine
    "57540" => 0.3, # NAD+
    "58885" => 0.3, # udp glucose
    "58297" => 0.3, # glutathione disulfide
    "58223" => 0.2, # UDP
    "57287" => 0.1, # coa
    "57692" => 0.01, # FAD
    "57783" => 0.01, # NADPH
    "57945" => 0.01, # NADH
    "58349" => 0.01, # nadph+
    "58210" => 0.01, # FMN
    "57986" => 0.001, # riboflavin
    "62501" => 0.001, # folate
    "597326" => 0.001, # pyridoxal 5 phosphate
    "60530" => 0.001, # heme o    
)

tot_sol = sum(values(solubles))
for (k, v) in solubles
    solubles[k] = v/tot_sol * molar_masses[k]
end

tot_sol = sum(values(solubles))
for (k, v) in solubles
    biomass[k] = get(biomass, k, 0) - soluble_pool * (v/tot_sol) * 1000 / molar_masses[k]
end

# required atp for growth
atp_req = 75.0

biomass["30616"] = biomass["30616"] - atp_req # atp
biomass["15377"] = -atp_req # water
biomass["43474"] = atp_req # pi
biomass["456216"] = atp_req # adp
biomass["15378"] = atp_req # h+

biomass

open(joinpath("data", "model", "biomass.json"), "w") do io
    JSON.print(io, biomass)
end
