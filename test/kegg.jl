@testset "KEGG utilities" begin

    stoich = VibrioNatriegens.parse_reaction_stoichiometry(
        "C00007 + 4 C00126 + 8 C00080 <=> 4 C00125 + 2 C00001 + 4 C00080",
    )
    @test stoich["C00080"] == -4 # -8 + 4
    @test stoich["C00126"] == -4
    @test stoich["C00001"] == 2


    r = VibrioNatriegens.get_kegg_reaction("R00081"; cache = true, force = true)

    m = VibrioNatriegens.get_kegg_compound("C00007"; cache = true, force = true)
    @test m.name == "Oxygen"
    @test m.formula == "O2"

    k = VibrioNatriegens.get_kegg_orthology("K01007"; cache = true, force = true)
    @test k.symbol == "pps, ppsA"
    @test first(k.reactions) == "R00199  ATP:pyruvate,water phosphotransferase"


    @test length(
        VibrioNatriegens.parse_kegg_definition(
            first(VibrioNatriegens.get_kegg_module("M00001")["DEFINITION"]),
        ),
    ) == 9
    @test length(
        VibrioNatriegens.parse_kegg_definition(
            first(VibrioNatriegens.get_kegg_module("M00151")["DEFINITION"]),
        ),
    ) == 1
    @test length(
        VibrioNatriegens.parse_kegg_definition(
            first(VibrioNatriegens.get_kegg_module("M00632")["DEFINITION"]),
        ),
    ) == 4
    @test length(
        VibrioNatriegens.parse_kegg_definition(
            first(VibrioNatriegens.get_kegg_module("M00005")["DEFINITION"]),
        ),
    ) == 1
    @test length(
        VibrioNatriegens.parse_kegg_definition(
            first(VibrioNatriegens.get_kegg_module("M00176")["DEFINITION"]),
        ),
    ) == 2


    # module completeness
    kegg_annos = CSV.File(
        joinpath("..", "data", "annotations", "kegg", "ko.txt");
        header = ["Protein", "KO", "Annotations"],
        delim = "\t",
    )
    all_kos = unique(
        String(ko) for x in [r for r in kegg_annos if !ismissing(r.KO)] for
        ko in split(x.KO, "/")
    )

    @test VibrioNatriegens.eval_module_completeness(
        VibrioNatriegens.parse_module_completeness("M00151"),
        all_kos,
    )
    @test VibrioNatriegens.eval_module_completeness(
        VibrioNatriegens.parse_module_completeness("M00001"),
        all_kos,
    )
    @test !VibrioNatriegens.eval_module_completeness(
        VibrioNatriegens.parse_module_completeness("M00580"),
        all_kos,
    )
    @test VibrioNatriegens.eval_module_completeness(
        VibrioNatriegens.parse_module_completeness("M00909"),
        all_kos,
    )
    @test !VibrioNatriegens.eval_module_completeness(
        VibrioNatriegens.parse_module_completeness("M00531"),
        all_kos,
    )
    @test VibrioNatriegens.eval_module_completeness(
        VibrioNatriegens.parse_module_completeness("M00016"),
        all_kos,
    )
    @test VibrioNatriegens.eval_module_completeness(
        VibrioNatriegens.parse_module_completeness("M00121"),
        all_kos,
    )

end
