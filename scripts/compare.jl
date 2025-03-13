using VibrioNatriegens
using AbstractFBCModels
using JSONFBCModels
using SBMLFBCModels
using COBREXA, Gurobi
using JSON
using ConstraintTrees
import ConstraintTrees as C
import AbstractFBCModels as A

model = VibrioNatriegens.build_model()
old = convert(A.CanonicalModel.Model, load_model("vnat_old.json"))

setdiff(A.reactions(old), A.reactions(model))

rids = filter(x -> isdigit(first(x)), A.reactions(old))
qts = get_quartets(rids)
lu = Dict(rr => r for (r, v) in qts for rr in v)
    
removed_reactions = []
stoichiometry_changed = []
bounds_changed = []
for rid in rids
    rid_new = lu[rid]
    rold = old.reactions[rid]
    lu[rid] in A.reactions(model) || (push!(removed_reactions, lu[rid]); continue)
    rnew =  model.reactions[lu[rid]]
    st_old = rold.stoichiometry
    st_new = rnew.stoichiometry
    length(st_old) == length(st_new) || @info(rid) 
    all(in.(collect(keys(st_old)), Ref(collect(keys(st_new))))) && all(in.(collect(keys(st_new)), Ref(collect(keys(st_old))))) || @info(rid)
    for (k, v) in st_old
        st_new[k] == v || push!(stoichiometry_changed, rid)
    end
    (rold.lower_bound != rnew.lower_bound || rold.upper_bound != rnew.upper_bound) && push!(bounds_changed, rid)
end

unique(removed_reactions)
