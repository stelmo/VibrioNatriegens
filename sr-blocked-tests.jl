using VibrioNatriegens
import AbstractFBCModels as A
using COBREXA, Gurobi
using JSON
using ConstraintTrees
import ConstraintTrees as C
using JSONFBCModels
using JSON

model = VibrioNatriegens.build_model()
flux_balance_analysis(model,optimizer=Gurobi.Optimizer)

#sol = flux_balance_analysis(model,optimizer=Gurobi.Optimizer)

###
ex_rids = filter(startswith("EX_"), A.reactions(model))
transporters = [rid for rid in A.reactions(model) if model.reactions[rid].transporter]
for rid in [ex_rids; transporters]
    model.reactions[rid].lower_bound = -1000
    model.reactions[rid].upper_bound = 1000
end
model.reactions["ATPM"].lower_bound = 0.0
###

# set metabolic.reactions to = -0.5 when needed rn, otherwise =-0.05
model.reactions["biomass"].stoichiometry["597326"] = -0.05 # pyridoxal
model.reactions["biomass"].stoichiometry["57287"] = -0.05 # coa
model.reactions["biomass"].stoichiometry["57920"] = -0.05 # norspermidine
model.reactions["biomass"].stoichiometry["57834"] = -0.05 # spermidine
model.reactions["biomass"].stoichiometry["18408"] = -0.05 # adenosylcobalamin
model.reactions["biomass"].stoichiometry["58437"] = -0.05 # deamido-nad
model.reactions["biomass"].stoichiometry["61386"] = -0.05 # d-alanine (UDP-N-acetyl-alpha-D-muramoyl-L-alanyl-gamma...)
model.reactions["biomass"].stoichiometry["67138"] = -0.05 # d-alanine (UDP-2-acetamido-2-deoxy-alpha-D-glucuronate)
model.reactions["biomass"].stoichiometry["65040"] = -0.05 # d-alanine (UDP-N-acetyl-alpha-D-galactosamine)
model.reactions["biomass"].stoichiometry["57451"] = -0.05 # 7,8-dihydrofolate
model.reactions["biomass"].stoichiometry["85130"] = -0.05 # 7,8-dihydroxanthopterin
model.reactions["biomass"].stoichiometry["62501"] = -0.05 # folate
model.reactions["biomass"].stoichiometry["60530"] = -0.5 # ferroheme o
model.reactions["biomass"].stoichiometry["58351"] = -0.5 # sirohydrochlorin
model.reactions["biomass"].stoichiometry["62414"] = -0.5 # biotinyl-5'-AMP, biotin pathway
model.reactions["biomass"].stoichiometry["64479"] = -0.5 # O-(pantetheine-4'-phosphoryl)-L-serine residue, biotin pathway
model.reactions["biomass"].stoichiometry["17790"] = -0.5 # methanol, biotin pathway
model.reactions["biomass"].stoichiometry["16490"] = -0.5 # S-adenosyl-4-methylsulfanyl-2-oxobutanoate, biotin pathway
model.reactions["biomass"].stoichiometry["59338"] = -0.5 # methylamine

model.reactions["EX_17154"].lower_bound = -1000.0 # set lb to -1000 to actually get nicotinamidase in
model.reactions["EX_58537"].lower_bound = -1000.0 # set lb to -1000 to actually get cob(II)yrinate a,c diamide in
model.reactions["EX_57761"].lower_bound = -1000.0 # set lb to -1000 to actually get pyridoxamine in
model.reactions["EX_16709"].lower_bound = -1000.0 # set lb to -1000 to actually get pyridoxine in
model.reactions["EX_2181"].lower_bound = -1000.0  # set lb to -1000 to actually get l-fucose in
model.reactions["EX_60344"].lower_bound = -1000.0 # set lb to -1000 to actually get heme in
model.reactions["EX_57433"].lower_bound = -1000.0 # set lb to -1000 to actually get sarcosine in
model.reactions["EX_57499"].lower_bound = -1000.0 # set lb to -1000 to actually get oxoadipate in
model.reactions["EX_58033"].lower_bound = -1000.0 # set lb to -1000 to actually get 2-phosphoglycolate in
model.reactions["EX_17742"].lower_bound = -1000.0 # set lb to -1000 to actually get 4-hydroxy-2-oxoglutarate in
model.reactions["EX_18391"].lower_bound = -1000.0 # set lb to -1000 to actually get d-gluconate in
model.reactions["EX_13941"].lower_bound = -1000.0 # set lb to -1000 to actually get carbamate in
model.reactions["EX_82735"].lower_bound = -1000.0 # set lb to -1000 to actually get pimeloyl-pantetheine-4-phosphorylserine methyl ester residue in, biotin pathway
model.reactions["EX_57360"].lower_bound = -1000.0 # set lb to -1000 to actually get 6-carboxyhexanoyl-CoA in, biotin pathway
model.reactions["EX_15354"].lower_bound = -1000.0 # set lb to -1000 to actually get choline in

model.reactions["sink_reaction1"] = VibrioNatriegens.Reaction( # need this
   stoichiometry = Dict("59338" => -1.0),
   lower_bound = 0.0,
   upper_bound = 1000.0,
)

model.reactions["sink_reaction2"] = VibrioNatriegens.Reaction( # need this
   stoichiometry = Dict("59338" => -1.0),
   lower_bound = 0.0,
   upper_bound = 1000.0,
)

model.reactions["sink_reaction3"] = VibrioNatriegens.Reaction( # need this
   stoichiometry = Dict("64479" => -1.0),
   lower_bound = 0.0,
   upper_bound = 1000.0,
)

#model.reactions["sink_reaction4"] = VibrioNatriegens.Reaction( # need this
#   stoichiometry = Dict("43474" => -1.0),
#   lower_bound = 0.0,
#   upper_bound = 1000.0,
#)

rid = "20217" # rhea id from blocked.csv 
rid = "sink_reaction1" # either rid with acutal id or the sink
ct = flux_balance_constraints(model)
sol1 = optimized_values(
    ct,
    objective = ct.fluxes[rid].value,
    sense = Minimal,
    optimizer = Gurobi.Optimizer,
)
sol2 = optimized_values(
    ct,
    objective = ct.fluxes[rid].value,
    sense = Maximal,
    optimizer = Gurobi.Optimizer,
)
lb = sol1.fluxes[rid]
ub = sol2.fluxes[rid]

 ct *= :parsimony^C.Constraint(
    sum(C.squared(x.value) for x in values(ct.fluxes)),
    nothing
)

ubct = deepcopy(ct)
ubct *= :ub^C.Constraint(
    ubct.fluxes[Symbol(rid)].value,
    C.EqualTo(ub),
)
sol2 = optimized_values(
    ubct,
    objective = ct.parsimony.value,
    sense = Minimal,
    optimizer = Gurobi.Optimizer,
)

open("vnat_fluxes.json", "w") do io
    JSON.print(io, sol2.fluxes)
end

C.pretty(
    C.ifilter_leaves(ub.fluxes) do ix, x
        abs(x) > 1e-6 && begin
            mets = [
                A.metabolite_name(model, k) for
                k in keys(A.reaction_stoichiometry(model, string(last(ix))))
            ]
            any(in.(mets, Ref(["pyridoxine"])))
        end
    end;
    format_label = x -> A.reaction_name(model, string(last(x))),
)