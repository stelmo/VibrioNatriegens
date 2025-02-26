using eQuilibrator, CSV, DataFrames, DataFramesMeta, RheaReactions, Measurements, Unitful

df = DataFrame(CSV.File(joinpath("data", "model", "metabolic_reactions.csv")))

eq = eQuilibrator.Equilibrator() # slow, only do if necessary

ds = Union{Missing,Float64}[]
pdgs = Union{Missing,Float64}[]

problem_rids = Int64[]

for r in eachrow(df)

    rid = r.RHEA_ID
    println(rid)

    coeff_mets = get_reaction_metabolites(rid)

    stoichiometry = Dict(string(v.accession) => s for (s, v) in coeff_mets)

    substrates = [string(abs(Int(v))) * " chebi:" * k for (k, v) in stoichiometry if v < 0]
    products = [string(abs(Int(v))) * " chebi:" * k for (k, v) in stoichiometry if v > 0]
    rxn_string = join(substrates, " + ") * " = " * join(products, " + ")

    d = missing
    pdg = missing
    try
        _d = ln_reversibility_index(eq, rxn_string; skip_unbalanced = true)
        _pdg = physiological_dg_prime(eq, rxn_string; skip_unbalanced = true)
        if !isnothing(_pdg)
            if Measurements.uncertainty(_pdg) > 0.5 * Measurements.value(_pdg) # error too big to be useful
                d = missing
                pdg = missing
            else
                d = Measurements.value(_d)
                pdg = Measurements.value(ustrip(u"kJ/mol", _pdg))
            end
        end
    catch

    end

    ismissing(d) && push!(problem_rids, rid)

    push!(ds, d)
    push!(pdgs, pdg)
end

unique(problem_rids)

# overwrite
df.DeltaG = pdgs
df.RevIndex = ds
df
CSV.write(joinpath("data", "model", "metabolic_reactions.csv"), df)
