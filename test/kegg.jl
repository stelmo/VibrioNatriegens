@testset "KEGG utilities" begin
    
    stoich = VibrioNatriegens.parse_reaction_stoichiometry("C00007 + 4 C00126 + 8 C00080 <=> 4 C00125 + 2 C00001 + 4 C00080")
    @test stoich["C00080"] == -4 # -8 + 4
    @test stoich["C00126"] == -4
    @test stoich["C00001"] == 2
    

    r = VibrioNatriegens.get_kegg_reaction("R00081"; cache=true, force=true)

    m = VibrioNatriegens.get_kegg_compound("C00007"; cache=true, force=true)
    @test m.name == "Oxygen"
    @test m.formula == "O2"

    k = VibrioNatriegens.get_kegg_orthology("K01007"; cache=true, force=true)
    @test k.symbol == "pps, ppsA"
    @test first(k.reactions) == "R00199  ATP:pyruvate,water phosphotransferase"

end
