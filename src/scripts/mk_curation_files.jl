
using CSV, DataFrames, DataFramesMeta, HTTP, VibrioNatriegens
using XLSX
using FASTX

prot = DataFrame(Protein = String[], Length = Int[])
FASTAReader(open(joinpath("data", "genome", "protein.faa"))) do reader
    for record in reader
        push!(prot, (identifier(record), length(sequence(record))))
    end
end

anno = DataFrame(CSV.File(joinpath("data", "curation", "all_annotations.csv")))
@rsubset!(anno, !ismissing(:KEGG_Reaction_Definition))
@rsubset!(anno, !ismissing(:EN_KO), !ismissing(:EN_Reaction))
anno = leftjoin(anno, prot, on = :Protein)

@select!(
    anno,
    :EN_Reaction,
    :Protein,
    :EN_EC,
    :EN_KO,
    :KEGG_KO,
    :HAMAP_Subunit,
    :Uniprot_Subunit,
    :EN_Symbol,
    :KEGG_Description,
    :RefSeq_Description,
    :KEGG_Reaction_Definition,
    :Length,
)

function guess_stoichiometry(s)
    ismissing(s) && return missing
    mer_map = Dict(
        "monomer" => 1,
        "homodimer" => 2,
        "homotetramer" => 4,
        "homotrimer" => 3,
        "homohexamer" => 6,
        "homopentamer" => 5,
        "homodecamer" => 10,
        "homooctamer" => 8,
        "homoheptamer" => 7,
        "homododecamer" => 12,
        "homomonomer" => 1,
    )
    ks = collect(keys(mer_map))
    lower_s = lowercase(s)
    for k in ks
        occursin(k, lower_s) && return mer_map[k]
    end
    return missing
end

function guess_stoichiometry(s1, s2)
    ss1 = guess_stoichiometry(s1)
    ss2 = guess_stoichiometry(s2)
    !ismissing(ss1) && return ss1
    !ismissing(ss2) && return ss2
    return 1 # default    
end

@rtransform!(
    anno,
    :GuessedStoichiometry = guess_stoichiometry(:HAMAP_Subunit, :Uniprot_Subunit)
)

all_kos = string.(Set(filter(!ismissing, [anno.EN_KO; anno.KEGG_KO])))
all_modules = VibrioNatriegens.list_kegg_modules()
complete_modules = String[]
for m in keys(all_modules)
    VibrioNatriegens.eval_module_completeness(
        VibrioNatriegens.parse_module_completeness(m),
        all_kos,
    ) && push!(complete_modules, m)
end

regex = r"R{1}(?:\d){5}"

# For complete modules only, find all their reactions with an annotated reaction

for modname in complete_modules
    println(modname)
    m = VibrioNatriegens.get_kegg_module(modname)
    zyme = []
    if length(m["REACTION"]) == 1 # these are isozymes, so all subunits are required
        rid = string(first(split(first(m["REACTION"]), "  ")))
        for orth in m["ORTHOLOGY"]
            kos, descr = string.(split(orth, "  "))
            for koterm in string.(split(kos, ','))
                ts = string.(split(string(koterm), '+'))
                all(in.(ts, Ref(all_kos))) && push!(zyme, ([rid], ts))
            end
        end
    else # these are multi-reaction modules
        for orth in m["ORTHOLOGY"]
            kos, descr = split(orth, " ", limit = 2)
            for koterm in split(kos, ',')
                ts = string.(split(string(koterm), '+'))
                if all(in.(ts, Ref(all_kos)))
                    rs = [string(x.match) for x in eachmatch(regex, descr)]
                    push!(zyme, (rs, ts))
                end
            end
        end
    end
    zyme

    XLSX.openxlsx(
        joinpath("data", "curation", "uncurated", "$modname.xlsx"),
        mode = "w",
    ) do xf
        n = length(zyme)
        for k = 1:n-1
            XLSX.addsheet!(xf)
        end
        sheet_names = String[]
        for (k, (rs, kos)) in enumerate(zyme)
            d = @rsubset(anno, any(skipmissing(in.(:EN_KO, Ref(kos)))))
            XLSX.writetable!(xf[k], d)
            sn = join(rs, "_")
            if sn in sheet_names
                kk = 'a'
                while true
                    if sn * kk in sheet_names
                        kk += 1
                    else
                        break
                    end
                end
                sn = sn * kk
            end

            push!(sheet_names, sn)
            XLSX.rename!(xf[k], sn)
        end
    end

end
