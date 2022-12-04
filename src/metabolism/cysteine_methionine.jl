module Cys_Met

using ..ModuleTools
using ..Utils

using COBREXA
using COBREXA.Types

subsystem = "Cysteine and Methionine Metabolism"

# (rhea_id = , name = "", isozymes = [[(),],], subsystem = subsystem,), #
 
rxns = (
    (rhea_id = 19169, name = "L-serine dehydratase", isozymes = [[(1, "A0A1B1ED13"),],[(1, "A0A1B1EHL2"),],[(1, "A0A1B1EJG8"),],], subsystem = subsystem,), #4.3.1.17

)

#=
These ECs were present but not included in rxns:


=#

end # module
