using eQuilibrator, CSV, DataFrames, DataFramesMeta, RheaReactions, Measurements, Unitful

df = DataFrame(CSV.File(joinpath("data", "model", "metabolic_reactions.csv")))

eq = eQuilibrator.Equilibrator() # slow, only do if necessary

ds = Union{Missing, Float64}[]
pdgs = Union{Missing, Float64}[]

problem_rids = Int64[]

for r in eachrow(df)
        
    rid = r.RHEA_ID
    println(rid)

    coeff_mets = get_reaction_metabolites(rid)

    stoichiometry = Dict(
        string(v.accession) => s
        for (s, v) in coeff_mets
    )

    substrates = [string(abs(Int(v))) * " chebi:" * k for (k, v) in stoichiometry if v < 0]
    products = [string(abs(Int(v))) * " chebi:" * k for (k, v) in stoichiometry if v > 0]
    rxn_string = join(substrates, " + ") * " = " * join(products, " + ")
    
    d = missing
    pdg = missing
    try 
        _d = ln_reversibility_index(eq, rxn_string; skip_unbalanced = true)
        _pdg =  physiological_dg_prime(eq, rxn_string; skip_unbalanced = true)
        d = isnothing(_d) ? missing : Measurements.value(_d)
        pdg = isnothing(_pdg) ? missing : Measurements.value(ustrip(u"kJ/mol", _pdg))
    catch
        
    end
    
    ismissing(d) && push!(problem_rids, rid)

    push!(ds, d)
    push!(pdgs, pdg)
end

# overwrite
df.DeltaG = pdgs
df.RevIndex = ds
df
CSV.write(joinpath("data", "model", "metabolic_reactions.csv"), df)
