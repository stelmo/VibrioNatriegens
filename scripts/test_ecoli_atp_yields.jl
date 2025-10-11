using AbstractFBCModels
import AbstractFBCModels as A
using Gurobi, COBREXA, ConstraintTrees
import ConstraintTrees as C
using JSONFBCModels

model = convert(A.CanonicalModel.Model, load_model("iML1515.json"))

model.reactions["BIOMASS_Ec_iML1515_core_75p37M"].objective_coefficient = 0.0
model.reactions["ATPM"].objective_coefficient = 1.0
model.reactions["ATPM"].lower_bound = 0.0
model.reactions["EX_glc__D_e"].lower_bound = 0.0 # glucose

csources = [ # (exchange rid, atp/substrate)
    (:EX_glc__D_e,) # glucose
    (:EX_ala__L_e, -) # alanine
    (:EX_succ_e,) # succinate
    (:EX_glyc_e,) # glycerol
    (:EX_glu__L_e,) # glutamate
    (:EX_rib__D_e,) # ribose
    (:EX_ac_e,) # acetate
    (:EX_gam_e,) # glucosamine
    (:EX_acgam_e,) # n-acetyl-d-glucosamine
    (:EX_glcn_e,) # gluconate
    (:EX_mal__L_e,) # malate
    (:EX_fum_e,) # fumarate
    (:EX_arab__L_e,) # L-arabinose
]
ct = flux_balance_constraints(model)
for (rid,) in csources
    ct.fluxes[rid].bound = C.Between(-10.0, 1000)
    sol = optimized_values(
        ct,
        optimizer = HiGHS.Optimizer,
        sense = Maximal,
        objective = ct.objective.value,
    )
    @info("ATP: $rid: $(sol.objective / sol.fluxes[rid])")
    ct.fluxes[rid].bound = C.Between(0.0, 1000)
end

# anaerobic
ct.fluxes.EX_glc__D_e.bound = C.Between(-10.0, 1000)
ct.fluxes.EX_o2_e.bound = C.EqualTo(0.0)
sol = optimized_values(
    ct,
    optimizer = HiGHS.Optimizer,
    sense = Maximal,
    objective = ct.objective.value,
)
@info("ATP: anaerobic: $(sol.objective / sol.fluxes[:EX_glc__D_e])")
