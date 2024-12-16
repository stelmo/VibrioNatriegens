using COBREXA, Gurobi
import AbstractFBCModels as A
import ConstraintTrees as C
using JSONFBCModels

ec = A.load(joinpath("data", "ecoli", "iML1515.json"))
sol = flux_balance_analysis(ec, optimizer = Gurobi.Optimizer)
C.pretty(C.filter_leaves(x -> abs(x) > 1, sol.fluxes))


