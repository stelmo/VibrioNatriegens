using COBREXA, Tulip, Clarabel
using Test 

include("code//reconstruction.jl")
import .Reconstruction as rc

# all tests assume the existence of the model
model = rc.build_model()

# run tests to ensure model is correct
@testset "Reconstruction tests" begin 
    for file in readdir("test//subtests")
        include("test//substests//$file")
    end
end