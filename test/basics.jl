
@testset "Mass balanced" begin
    
    model = VibrioNatriegens.build_model()

    rids = filter(x -> !startswith(x, "EX_") && x != "biomass", A.reactions(model))
    unbal_rids = String[]
    for rid in rids
        s = A.reaction_stoichiometry(model, rid)
        m = Dict()
        for (k, v) in s
            for (kk, vv) in A.metabolite_formula(model, k)
                m[kk] = get(m, kk, 0) + vv * v 
            end
        end
        all(values(m) .== 0) || push!(unbal_rids, rid)    
    end

    @test isempty(unbal_rids)
end

@testset "Charge balanced" begin
    
    model = VibrioNatriegens.build_model()

    rids = filter(x -> !startswith(x, "EX_") && x != "biomass", A.reactions(model))
    unbal_rids = String[]
    for rid in rids
        s = A.reaction_stoichiometry(model, rid)
        m = 0
        for (k, v) in s
            m += v * A.metabolite_charge(model, k) 
        end
        m == 0 || push!(unbal_rids, rid)    
    end

    @test isempty(unbal_rids)
end

@testset "ATP production" begin

    model = VibrioNatriegens.build_model()

    model.reactions["biomass"].objective_coefficient = 0.0
    model.reactions["ATPM"].objective_coefficient = 1.0
    model.reactions["ATPM"].lower_bound = 0.0
    

    sol = flux_balance_analysis(model, optimizer=Gurobi.Optimizer)
    llsol = loopless_flux_balance_analysis(model, optimizer=Gurobi.Optimizer)
    @test isapprox(sol.objective, llsol.objective)
    
    for rid in filter(x -> startswith(x, "EX_"), A.reactions(model))
        model.reactions[rid].lower_bound = 0.0
        model.reactions[rid].upper_bound = 0.0
    end

    sol = flux_balance_analysis(model, optimizer=Gurobi.Optimizer)
    @test isapprox(sol.objective, 0.0)
     
end
