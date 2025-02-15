

function is_mass_balanced(model, rid)

    stoic = A.reaction_stoichiometry(model, rid)
    atoms = Dict()
    for (k, v) in stoic
        for (kk, vv) in A.metabolite_formula(model, k)
            atoms[kk] = get(atoms, kk, 0) + v * vv
        end
    end
    
    all(values(atoms) .== 0)
end
