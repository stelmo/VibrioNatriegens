@testset "Overview" begin 
    sol = flux_balance_analysis_dict(model, Tulip.Optimizer)
    
    # can grow in defined media
    @test sol["biomass"] ≈ 0.5

end