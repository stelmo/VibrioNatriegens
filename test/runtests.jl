using COBREXA
using Test 

include("code//reconstruction.jl")
import .Reconstruction as rc

model = rc.build_model()

# run tests to ensure model is correct
